#include <algorithm>
#include <chrono>
#include <cmath>
#include <execution>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

const double mu = 1.0;
const double alpha = 2.0;
const double bar_n = 0.9;
const double p = 1.1;
const double c = 1.0E-9;
const double beta = std::log(10);

const int bar_width = 40;
const double delay = 0.5;    // refresh delay in seconds

void update_bar(const double progress, const std::string& status = "")
{
    const int pos = static_cast<int>(progress * bar_width);

    std::cout << '[';
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos)        std::cout << '=';
        else if (i == pos)  std::cout << '>';
        else                std::cout << ' ';
    }
    std::cout << "] " << std::setw(3)
        << static_cast<int>(progress * 100) << "%";
    if (status != "")
        std::cout << " -- " << status;
    std::cout << '\r';
    std::cout.flush();
}

void close_bar(const std::string& status = "")
{
    update_bar(1, status);
    std::cout << '\n';
}

std::string to_string(const double x, int precision = 2)
{
    std::ostringstream ss;

    ss << std::setprecision(precision) << x;
    return ss.str();
}

namespace rd
{
    inline std::mt19937 init()
    {
        std::random_device rd;
        std::seed_seq ss{ rd(), rd(), rd(), rd() };

        return std::mt19937{ ss };
    }

    inline std::mt19937 mt = init();
    inline std::uniform_real_distribution<double> d{ 0, 1 };

    inline double logrand()
    {
        return std::log(d(mt));
    }

    inline double mag()
    {
        return -1 / beta * logrand();
    }
}

struct Point
{
    double          t;
    double          m;
    std::size_t     parent;
};

using Sequence = std::vector<Point>;

Sequence etas(const double tmax = 100.0, const bool verbose = true)
{
    Sequence seq;
    double tc = 0;
    double m_max = 0;

    if (verbose)
        std::cout << "Generating background earthquakes...\n";

    while (tc < tmax) {
        const double dt = -1 / mu * rd::logrand();

        tc += dt;
        if (tc < tmax) {
            const double m = rd::mag();

            m_max = std::max(m_max, m);
            seq.push_back({ tc, m, 0 });
        }
    }

    const double a = bar_n * (p - 1) * (beta - alpha)
        / (beta * std::pow(c, 1 - p));

    if (a > 0) {
        std::size_t nc = 0;
        auto start = std::chrono::high_resolution_clock::now();

        if (verbose) {
            std::cout << "Generating aftershocks...\n";
            update_bar(seq[nc].t / tmax, to_string(m_max));
        }

        while (true) {
            tc = 0;

            if (verbose) {
                const auto end = std::chrono::high_resolution_clock::now();
                const std::chrono::duration<double> diff = end - start;

                if (diff.count() > delay) {
                    update_bar(seq[nc].t / tmax, to_string(m_max));
                    start = end;
                }
            }

            while (true) {
                const double tmp = std::pow(tc + c, 1 - p) + (p - 1)
                    / (a * std::exp(alpha * seq[nc].m))
                    * rd::logrand();

                if (tmp > 0) {
                    const double dt = std::pow(tmp, 1 / (1 - p)) - tc - c;
                    tc += dt;
                    const double tc_nc = tc + seq[nc].t;

                    if (tc_nc < tmax) {
                        const double m = rd::mag();

                        m_max = std::max(m_max, m);
                        seq.push_back({ tc_nc, m, nc + 1 });
                    }
                    else {
                        break;
                    }
                }
                else {
                    break;
                }
            }

            auto sort = [](const Point& p1, const Point& p2)
            {
                return p1.t < p2.t;
            };

            std::sort(std::execution::par, seq.begin(), seq.end(), sort);

            ++nc;
            if (nc >= seq.size())
                break;
        }
    }

    if (verbose)
        close_bar(to_string(m_max));

    return seq;
}

void write_to_file(const Sequence& seq, const std::string& filename,
    bool verbose = true)
{
    std::ofstream file{ filename };
    std::size_t id = 0;

    file << "ID,TIME,MAG,PARENT\n";
    for (const Point& e : seq) {
        file << ++id << ',' << e.t << ',' << e.m << ',' << e.parent << '\n';
    }

    if (verbose)
        std::cout << seq.size() << " events written to file `"
        << filename << "`.\n";
}

void generate_seqs(int num_seqs, std::size_t max_len, double tmax = 100.0,
    const std::string& dirname = "data", bool verbose = true)
{
    if (!std::filesystem::create_directory(dirname)) {
        std::cout << "Could not create directory, exiting...\n";
        return;
    }

    if (verbose)
        std::cout << "Generating sequences...\n";

    for (int i = 1; i <= num_seqs; ++i) {
        do {
            const Sequence seq = etas(tmax, false);

            if (seq.size() <= max_len) {
                const std::string filename = dirname + "/data"
                    + std::to_string(i) + ".csv";

                write_to_file(seq, filename, false);
                break;
            }
        } while (true);

        if (verbose)
            update_bar(static_cast<double>(i) / num_seqs);
    }

    if (verbose)
        close_bar();
}

int main()
{
    const std::string dirname = "data";
    double tmax = 100.0;
    generate_seqs(100, 300, tmax, dirname);
    
    const auto tmax50 = etas(50.0);
    write_to_file(tmax50, "tmax50.csv");
    
    const auto tmax5k = etas(5000.0);
    write_to_file(tmax5k, "tmax5k.csv");

    return 0;
}
