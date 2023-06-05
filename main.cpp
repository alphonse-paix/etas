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

#include <boost/program_options.hpp>
namespace po = boost::program_options;

namespace args {
    double mu = 1;
    double alpha = 2;
    double bar_n = 0.9;
    double p = 1.1;
    double c = 1e-09;
    double beta = std::log(10);
    double tmax = 1e+03;
    std::string filename = "data.csv";
    std::string dirname = "data";
}

const int bar_width = 40;
const double delay = 0.5;    // refresh delay in seconds

void update_bar(const double progress, const std::string& status = "") {
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

void close_bar(const std::string& status = "") {
    update_bar(1, status);
    std::cout << '\n';
}

std::string to_string(const double x, int precision = 2) {
    std::ostringstream ss;

    ss << std::setprecision(precision) << x;
    return ss.str();
}

namespace rd {
    inline std::mt19937 init() {
        std::random_device rd;
        std::seed_seq ss{ rd(), rd(), rd(), rd() };

        return std::mt19937{ ss };
    }

    inline std::mt19937 mt = init();
    inline std::uniform_real_distribution<double> d{ 0, 1 };

    inline double logrand() {
        return std::log(d(mt));
    }

    inline double mag(double beta) {
        return -1 / beta * logrand();
    }
}

struct Point {
    double          t;
    double          m;
    std::size_t     parent;
};

using Sequence = std::vector<Point>;

Sequence etas(const bool verbose = true) {
    Sequence seq;
    double tc = 0;
    double m_max = 0;

    if (verbose)
        std::cout << "Generating background earthquakes...\n";

    while (tc < args::tmax) {
        const double dt = -1 / args::mu * rd::logrand();

        tc += dt;
        if (tc < args::tmax) {
            const double m = rd::mag(args::beta);

            m_max = std::max(m_max, m);
            seq.push_back({ tc, m, 0 });
        }
    }

    const double a = args::bar_n * (args::p - 1) * (args::beta - args::alpha)
        / (args::beta * std::pow(args::c, 1 - args::p));

    if (a > 0) {
        std::size_t nc = 0;
        auto start = std::chrono::high_resolution_clock::now();

        if (verbose) {
            std::cout << "Generating aftershocks...\n";
            update_bar(seq[nc].t / args::tmax, to_string(m_max));
        }

        while (true) {
            tc = 0;

            if (verbose) {
                const auto end = std::chrono::high_resolution_clock::now();
                const std::chrono::duration<double> diff = end - start;

                if (diff.count() > delay) {
                    update_bar(seq[nc].t / args::tmax, to_string(m_max));
                    start = end;
                }
            }

            while (true) {
                const double tmp = std::pow(tc + args::c, 1 - args::p) + (args::p - 1)
                    / (a * std::exp(args::alpha * seq[nc].m))
                    * rd::logrand();

                if (tmp > 0) {
                    const double dt = std::pow(tmp, 1 / (1 - args::p)) - tc - args::c;
                    tc += dt;
                    const double tc_nc = tc + seq[nc].t;

                    if (tc_nc < args::tmax) {
                        const double m = rd::mag(args::beta);

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
    bool verbose = true) {
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

void generate_seqs(int num_seqs, std::size_t max_len,
    const std::string& dirname = "data", bool verbose = true) {
    if (!std::filesystem::create_directory(dirname)) {
        std::cout << "Could not create directory, exiting...\n";
        return;
    }

    if (verbose)
        std::cout << "Generating sequences...\n";

    for (int i = 1; i <= num_seqs; ++i) {
        do {
            const Sequence seq = etas(false);

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

void print_args() {
    std::cout << "Executing program with following parameters:\n";
    std::cout << "tmax\t\t" << args::tmax << '\n';
    std::cout << "mu\t\t" << args::mu << '\n';
    std::cout << "alpha\t\t" << args::alpha << '\n';
    std::cout << "bar_n\t\t" << args::bar_n << '\n';
    std::cout << "p\t\t" << args::p << '\n';
    std::cout << "c\t\t" << args::c << '\n';
    std::cout << "beta\t\t" << args::beta << '\n';
    std::cout << "filename\t" << args::filename << '\n';
    std::cout << "dirname\t\t" << args::dirname << '\n';
}

int main(int argc, char* argv[]) {
    // arguments parsing
    // -----------------
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "ETAS algorithm implementation")
        ("tmax", po::value<double>(), "")
        ("mu", po::value<double>(), "")
        ("alpha", po::value<double>(), "")
        ("bar_n", po::value<double>(), "")
        ("p", po::value<double>(), "")
        ("c", po::value<double>(), "")
        ("beta", po::value<double>(), "")
        ("filename", po::value<std::string>(), "")
        ("dirname", po::value<std::string>(), "")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 0;
    }

    if (vm.count("tmax")) {
        std::cout << "Setting tmax...\n";
        args::tmax = vm["tmax"].as<double>();
    }

    if (vm.count("mu")) {
        std::cout << "Setting mu...\n";
        args::mu = vm["mu"].as<double>();
    }

    if (vm.count("alpha")) {
        std::cout << "Setting alpha...\n";
        args::alpha = vm["alpha"].as<double>();
    }

    if (vm.count("bar_n")) {
        std::cout << "Setting bar_n...\n";
        args::bar_n = vm["bar_n"].as<double>();
    }

    if (vm.count("p")) {
        std::cout << "Setting p...\n";
        args::p = vm["p"].as<double>();
    }

    if (vm.count("c")) {
        std::cout << "Setting c...\n";
        args::c = vm["c"].as<double>();
    }

    if (vm.count("beta")) {
        std::cout << "Setting beta...\n";
        args::beta = vm["beta"].as<double>();
    }

    if (vm.count("filename")) {
        std::cout << "Setting filename...\n";
        args::filename = vm["filename"].as<std::string>();
    }

    if (vm.count("dirname")) {
        std::cout << "Setting dirname...\n";
        args::dirname = vm["dirname"].as<std::string>();
    }

    // print ETAS simulation arguments
    // -------------------------------
    print_args();

    // simulation
    // ----------    
    const auto tmax3k = etas();
    write_to_file(tmax3k, args::filename);

    return 0;
}
