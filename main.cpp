#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
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

const int bar_width = 40;
const double delay = 0.5;   // progress bar update delay in s

struct Config {
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

Sequence etas(const Config& cfg) {
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
            update_bar(seq[nc].t / cfg.tmax, to_string(m_max));
        }

        while (true) {
            tc = 0;

            if (verbose) {
                const auto end = std::chrono::high_resolution_clock::now();
                const std::chrono::duration<double> diff = end - start;

                if (diff.count() > delay) {
                    update_bar(seq[nc].t / cfg.tmax, to_string(m_max));
                    start = end;
                }
            }

            while (true) {
                const double tmp = std::pow(tc + cfg.c, 1 - cfg.p) + (cfg.p - 1)
                    / (a * std::exp(cfg.alpha * seq[nc].m))
                    * rd::logrand();

                if (tmp > 0) {
                    const double dt = std::pow(tmp, 1 / (1 - cfg.p)) - tc - cfg.c;
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

            auto sort = [](const Point& p1, const Point& p2) {
                return p1.t < p2.t;
            };

            std::sort(std::execution::par, seq.begin(), seq.end(), sort);

            ++nc;
            if (nc >= seq.size())
                break;

            if (seq.size() > cfg.max_len && cfg.generate_seqs)
                return seq;     // early return, will not be saved
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

void generate_seqs(const Config& cfg) {
    if (!std::filesystem::create_directory(cfg.dirname)) {
        std::cout << "Could not create directory, exiting...\n";
        return;
    }

    if (cfg.verbose)
        std::cout << "Generating sequences...\n";

    for (int i = 1; i <= cfg.num_seqs; ++i) {
        do {
            const Sequence seq = etas(cfg);

            if (seq.size() <= cfg.max_len) {
                const std::string filename = cfg.dirname + '/' + cfg.filename
                    + std::to_string(i) + ".csv";

                write_to_file(seq, filename, false);
                break;
            }
        } while (true);

        if (cfg.verbose)
            update_bar(static_cast<double>(i) / cfg.num_seqs);
    }

    if (cfg.verbose)
        close_bar();
}

Config parse_arguments(int argc, char* argv[]) {
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "")
        ("tmax", po::value<double>()->default_value(100), "")
        ("mu", po::value<double>()->default_value(1), "")
        ("alpha", po::value<double>()->default_value(2), "")
        ("bar_n", po::value<double>()->default_value(0.9), "")
        ("p", po::value<double>()->default_value(1.1), "")
        ("c", po::value<double>()->default_value(1e-09), "")
        ("beta", po::value<double>()->default_value(std::log(10)), "")

        ("generate_seqs", po::bool_switch()->default_value(false), "")
        ("num_seqs", po::value<int>()->default_value(100), "")
        ("max_len", po::value<int>()->default_value(300), "")

        ("filename", po::value<std::string>()->default_value("data.csv"), "")
        ("dirname", po::value<std::string>()->default_value("data"), "")

        ("verbose", po::bool_switch()->default_value(false), "")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << '\n';
        std::exit(0);
    }

    return Config{
        vm["tmax"].as<double>(),
        vm["mu"].as<double>(),
        vm["alpha"].as<double>(),
        vm["bar_n"].as<double>(),
        vm["p"].as<double>(),
        vm["c"].as<double>(),
        vm["beta"].as<double>(),

        vm["generate_seqs"].as<bool>(),
        vm["num_seqs"].as<int>(),
        vm["max_len"].as<int>(),

        vm["filename"].as<std::string>(),
        vm["dirname"].as<std::string>(),

        vm["verbose"].as<bool>()
    };
}

void print_args(const Config& cfg) {
    std::cout << "Executing program with following parameters:\n";
    std::cout << std::boolalpha;

    std::cout << "tmax\t\t" << cfg.tmax << '\n';
    std::cout << "mu\t\t" << cfg.mu << '\n';
    std::cout << "alpha\t\t" << cfg.alpha << '\n';
    std::cout << "bar_n\t\t" << cfg.bar_n << '\n';
    std::cout << "p\t\t" << cfg.p << '\n';
    std::cout << "c\t\t" << cfg.c << '\n';
    std::cout << "beta\t\t" << cfg.beta << "\n\n";

    std::cout << "generate_seqs\t" << cfg.generate_seqs << '\n';
    std::cout << "num_seqs\t" << cfg.num_seqs << '\n';
    std::cout << "max_len\t\t" << cfg.max_len << "\n\n";

    std::cout << "filename\t" << cfg.filename << '\n';
    std::cout << "dirname\t\t" << cfg.dirname << "\n\n";

    std::cout << "verbose\t\t" << cfg.verbose << '\n';
}

int main(int argc, char* argv[]) {
    const Config cfg = parse_arguments(argc, argv);

    if (cfg.verbose)
        print_args(cfg);
    
    if (cfg.generate_seqs) {
        generate_seqs(cfg);
    }
    else {
        auto seq = etas(cfg);
        write_to_file(seq, cfg.filename, cfg.verbose);
    }

    return 0;
}
