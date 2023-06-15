#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <execution>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "argparse.hpp"

namespace pbar_settings
{
    constexpr int bar_width     = 40;
    constexpr double delay      = 0.5;
    constexpr char fill         = '=';
    constexpr char lead         = '>';
    constexpr char remainder    = ' ';
}

struct Config
{
    double tmax;
    double mu;
    double alpha;
    double bar_n;
    double p;
    double c;
    double beta;

    bool generate_seqs;
    int num_seqs;
    int max_len;

    std::string filename;
    std::string dirname;

    bool verbose;
};

template<typename T>
void update_bar(const double progress, const T& status)
{
    const int pos = static_cast<int>(progress * pbar_settings::bar_width);

    std::cout << '[';
    for (int i = 0; i < pbar_settings::bar_width; ++i) {
        if (i < pos)        std::cout << pbar_settings::fill;
        else if (i == pos)  std::cout << pbar_settings::lead;
        else                std::cout << pbar_settings::remainder;
    }
    std::cout << "] " << std::setw(3) << static_cast<int>(progress * 100)
                      << "% " << status;
    std::cout << '\r';
    std::cout.flush();
}

template<typename T>
void close_bar(const T& status)
{
    update_bar(1, status);
    std::cout << '\n';
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

    inline double mag(double beta)
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

Sequence etas(const Config& cfg)
{
    Sequence seq;
    double tc = 0;
    double m_max = 0;
    bool verbose = cfg.verbose && !cfg.generate_seqs;

    if (verbose)
        std::cout << "Generating background earthquakes...\n";

    while (tc < cfg.tmax) {
        const double dt = -1 / cfg.mu * rd::logrand();

        tc += dt;
        if (tc < cfg.tmax) {
            const double m = rd::mag(cfg.beta);

            m_max = std::max(m_max, m);
            seq.push_back({ tc, m, 0 });
        }
    }

    const double a = cfg.bar_n * (cfg.p - 1) * (cfg.beta - cfg.alpha)
        / (cfg.beta * std::pow(cfg.c, 1 - cfg.p));

    if (a > 0) {
        std::size_t nc = 0;
        auto start = std::chrono::high_resolution_clock::now();

        if (verbose) {
            std::cout << "Generating aftershocks...\n";
            update_bar(seq[nc].t / cfg.tmax, m_max);
        }

        while (true) {
            tc = 0;

            if (verbose) {
                const auto end = std::chrono::high_resolution_clock::now();
                const std::chrono::duration<double> diff = end - start;

                if (diff.count() > pbar_settings::delay) {
                    update_bar(seq[nc].t / cfg.tmax, m_max);
                    start = end;
                }
            }

            while (true) {
                const double tmp = std::pow(tc + cfg.c, 1 - cfg.p)
                    + (cfg.p - 1) / (a * std::exp(cfg.alpha * seq[nc].m))
                    * rd::logrand();

                if (tmp > 0) {
                    const double dt = std::pow(tmp, 1 / (1 - cfg.p))
                        - tc - cfg.c;
                    tc += dt;
                    const double tc_nc = tc + seq[nc].t;

                    if (tc_nc < cfg.tmax) {
                        const double m = rd::mag(cfg.beta);

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

            if (static_cast<int>(seq.size()) > cfg.max_len
                    && cfg.generate_seqs)
                return seq;
        }
    }

    if (verbose)
        close_bar(m_max);

    return seq;
}

void write_to_file(const Sequence& seq,
                   const std::string& filename,
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

void generate_seqs(const Config& cfg)
{
    if (!std::filesystem::create_directory(cfg.dirname)) {
        std::cout << "Could not create directory, exiting...\n";
        return;
    }

    if (cfg.verbose)
        std::cout << "Generating sequences...\n";

    for (int i = 1; i <= cfg.num_seqs; ++i) {
        do {
            const Sequence seq = etas(cfg);

            if (static_cast<int>(seq.size()) <= cfg.max_len) {
                const std::string filename = cfg.dirname + '/' + cfg.filename
                    + std::to_string(i) + ".csv";

                write_to_file(seq, filename, false);
                break;
            }
        } while (true);

        if (cfg.verbose)
            update_bar(static_cast<double>(i) / cfg.num_seqs, "");
    }

    if (cfg.verbose)
        close_bar("");
}

void print_args(const Config& cfg)
{
    std::cout << "Program configuration:\n";
    std::cout << std::boolalpha;

    std::cout << "  --tmax\t\t" << cfg.tmax << '\n';
    std::cout << "  --mu\t\t\t" << cfg.mu << '\n';
    std::cout << "  --alpha\t\t" << cfg.alpha << '\n';
    std::cout << "  --bar_n\t\t" << cfg.bar_n << '\n';
    std::cout << "  --p\t\t\t" << cfg.p << '\n';
    std::cout << "  --c\t\t\t" << cfg.c << '\n';
    std::cout << "  --beta\t\t" << cfg.beta << "\n\n";

    std::cout << "  --generate_seqs\t" << cfg.generate_seqs << '\n';
    std::cout << "  --num_seqs\t\t" << cfg.num_seqs << '\n';
    std::cout << "  --max_len\t\t" << cfg.max_len << "\n\n";

    std::cout << "  --filename\t\t" << cfg.filename << '\n';
    std::cout << "  --dirname\t\t" << cfg.dirname << "\n\n";

    std::cout << "  --verbose\t\t" << cfg.verbose << '\n';
}

Config parse_arguments(int argc, char* argv[])
{
    argparse::ArgumentParser program{ "ETAS" };

    program.add_argument("--tmax").default_value(100.0).scan<'g', double>();
    program.add_argument("--mu").default_value(1.0).scan<'g', double>();
    program.add_argument("--alpha").default_value(2.0).scan<'g', double>();
    program.add_argument("--bar_n").default_value(0.9).scan<'g', double>();
    program.add_argument("--p").default_value(1.1).scan<'g', double>();
    program.add_argument("--c").default_value(1e-09).scan<'g', double>();
    program.add_argument("--beta").default_value(std::log(10))
        .scan<'g', double>();

    program.add_argument("--generate_seqs").default_value(false)
        .implicit_value(true);
    program.add_argument("--num_seqs").default_value(100)
        .scan<'i', int>();
    program.add_argument("--max_len").default_value(300)
        .scan<'i', int>();

    program.add_argument("--filename")
        .default_value(std::string{ "data.csv" });
    program.add_argument("--dirname").default_value(std::string{ "data" });

    program.add_argument("--verbose").default_value(false)
        .implicit_value(true);

    program.add_argument("--print_config").default_value(false)
        .implicit_value(true);

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        std::cerr << err.what() << '\n';
        std::cerr << program;
        std::exit(1);
    }

    Config cfg{
        program.get<double>("--tmax"),
        program.get<double>("--mu"),
        program.get<double>("--alpha"),
        program.get<double>("--bar_n"),
        program.get<double>("--p"),
        program.get<double>("--c"),
        program.get<double>("--beta"),

        program.get<bool>("--generate_seqs"),
        program.get<int>("--num_seqs"),
        program.get<int>("--max_len"),

        program.get<std::string>("--filename"),
        program.get<std::string>("--dirname"),

        program.get<bool>("--verbose"),
    };

    if (program["--print_config"] == true)
        print_args(cfg);

    return cfg;
}

int main(int argc, char* argv[])
{
    const Config cfg = parse_arguments(argc, argv);

    if (cfg.generate_seqs)
        generate_seqs(cfg);
    else
        write_to_file(etas(cfg), cfg.filename, cfg.verbose);

    return 0;
}
