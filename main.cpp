// main.cpp
// ETAS algorithm
// argpare, indicators (github.com/p-ranav)

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

#include <argparse/argparse.hpp>
#include <indicators/progress_bar.hpp>
#include <indicators/cursor_control.hpp>

// forward declarations
// ------------------------------------------------------------------------------------------------

struct Config;
Config parse_arguments(int argc, char* argv[]);
void print_args(const Config& cfg);

/*
namespace rd
    double logrand()
    dobule mag(double beta)
*/

struct Event;
using Sequence = std::vector<Event>;
Sequence etas(const Config& cfg);
void generate_seqs(const Config& cfg);
void write_to_file(const Sequence& seq, const std::string& filename, bool verbose);

// program configuration related functions
// ------------------------------------------------------------------------------------------------

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

// parsing command line arguments with argparse from p-ranav
// ---------------------------------------------------------
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
        .default_value(std::string{ "data" });
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

// print arguments if print_config is passed
// -----------------------------------------
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

// progress bar
// ------------------------------------------------------------------------------------------------

constexpr double refresh_delay_in_s = 0.5;

auto get_progress_bar(const std::string_view& prefix_text) -> indicators::ProgressBar
{
    using namespace indicators;

    return ProgressBar{
        option::BarWidth{40},
        option::Start{"["},
        option::Fill{"="},
        option::Lead{">"},
        option::Remainder{" "},
        option::End{"]"},
        option::ForegroundColor{Color::white},
        option::ShowElapsedTime{true},
        option::PrefixText{prefix_text},
        option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
    };
}

// simulation
// ------------------------------------------------------------------------------------------------

// random functions required for the simulation
// --------------------------------------------
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

struct Event
{
    double          t;
    double          m;
    std::size_t     parent;
};

// main logic for the ETAS algorithm
// ---------------------------------
Sequence etas(const Config& cfg)
{
    using namespace indicators;

    Sequence seq;
    double tc = 0;
    double m_max = 0;
    bool verbose = cfg.verbose && !cfg.generate_seqs;
    auto progress_bar = get_progress_bar("Generating aftershocks...");

    if (verbose)
        std::cout << "Generating background earthquakes... \n";

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

    if (verbose) {
        progress_bar.set_progress(0);
        progress_bar.set_option(option::PostfixText{ std::to_string(m_max) });
    }

    if (a > 0) {
        std::size_t nc = 0;
        auto start = std::chrono::high_resolution_clock::now();

        while (true) {
            tc = 0;

            if (verbose) {
                const auto end = std::chrono::high_resolution_clock::now();
                const std::chrono::duration<double> diff = end - start;

                if (diff.count() > refresh_delay_in_s) {
                    progress_bar.set_progress(seq[nc].t / cfg.tmax * 100.0);
                    progress_bar.set_option(option::PostfixText{ std::to_string(m_max) });
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

            auto sort = [](const Event& p1, const Event& p2)
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
        progress_bar.set_progress(100);

    return seq;
}

// generate sequences and write each to its own file
// -------------------------------------------------
void generate_seqs(const Config& cfg)
{
    if (!std::filesystem::create_directory(cfg.dirname)) {
        std::cout << "Could not create directory, exiting...\n";
        return;
    }

    auto progress_bar = get_progress_bar("Generating sequences... ");
    progress_bar.set_progress(0);

    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 1; i <= cfg.num_seqs; ++i) {
        do {
            const Sequence seq = etas(cfg);

            if (static_cast<int>(seq.size()) <= cfg.max_len) {
                const std::string filename = cfg.dirname + '/' + cfg.filename
                    + std::to_string(i);

                write_to_file(seq, filename, false);
                break;
            }
        } while (true);

        const auto end = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> diff = end - start;

        if (cfg.verbose) {
            if (diff.count() > refresh_delay_in_s) {
                progress_bar.set_progress(static_cast<double>(i) / cfg.num_seqs * 100.0);
                start = end;
            }
        }
    }

    progress_bar.set_progress(100);
}

// take a sequence, a filename and write that sequence to the specified file with a CSV format
// -------------------------------------------------------------------------------------------
void write_to_file(const Sequence& seq,
    const std::string& filename,
    bool verbose = true)
{
    std::ofstream file{ filename + ".csv" };
    std::size_t id = 0;

    file << "ID,TIME,MAG,PARENT\n";
    for (const Event& e : seq) {
        file << ++id << ',' << e.t << ',' << e.m << ',' << e.parent << '\n';
    }

    if (verbose)
        std::cout << seq.size() << " events written to file `"
        << filename << ".csv`.\n";
}

// ------------------------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    const Config cfg = parse_arguments(argc, argv);

    indicators::show_console_cursor(false);

    if (cfg.generate_seqs)
        generate_seqs(cfg);
    else
        write_to_file(etas(cfg), cfg.filename, cfg.verbose);

    indicators::show_console_cursor(true);
    return 0;
}
