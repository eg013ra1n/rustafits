//! `rustafits` — FITS/XISF astronomical images from the terminal: convert to
//! auto-stretched JPEG/PNG (the default command), list headers, measure stars.

mod analyze;
mod batch;
mod convert;
mod info;
mod report;

use clap::{CommandFactory, Parser, Subcommand};

const EXAMPLES: &str = "\
Examples:
  rustafits light.fits                      # → light.jpg next to it
  rustafits m31.fits m31.png                # classic form: input, output
  rustafits *.fits -o previews/             # batch into a directory
  rustafits night1/ -r -o out/ -n           # a tree, mirrored; skip existing
  rustafits m31.xisf -o - | viewer -        # JPEG on stdout
  rustafits info night1/                    # exposure, filter, gain, temp per file
  rustafits analyze night1/ --sort fwhm     # stars, FWHM, HFR, eccentricity, SNR
  rustafits completions zsh > ~/.zfunc/_rustafits";

#[derive(Parser, Debug)]
#[command(
    name = "rustafits",
    version,
    about = "Convert FITS/XISF astronomical images to auto-stretched JPEG/PNG, list headers, measure stars",
    after_help = EXAMPLES,
    args_conflicts_with_subcommands = true,
    subcommand_negates_reqs = true
)]
struct Cli {
    #[command(subcommand)]
    command: Option<Command>,

    #[command(flatten)]
    convert: convert::ConvertArgs,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Convert to JPEG/PNG — the default, so `convert` may be omitted (use it
    /// for an input file literally named like a subcommand)
    Convert(convert::ConvertArgs),
    /// Header summary without reading pixels: size, type, exposure, filter, gain, temperature
    Info(info::InfoArgs),
    /// Star metrics per frame: stars, FWHM, HFR, eccentricity, SNR (library defaults, as in Athenaeum)
    Analyze(analyze::AnalyzeArgs),
    /// Print a shell completion script (bash, zsh, fish, elvish, powershell)
    Completions {
        #[arg(value_enum)]
        shell: clap_complete::Shell,
    },
}

fn main() {
    // Rust ignores SIGPIPE, so `rustafits info … | head` would panic on the
    // first write after `head` exits. Restore the default: exit quietly.
    #[cfg(unix)]
    // SAFETY: called first thing in main, before any other thread exists.
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }

    let cli = Cli::parse();
    let result = match cli.command {
        None => convert::run(cli.convert),
        Some(Command::Convert(args)) => convert::run(args),
        Some(Command::Info(args)) => info::run(args),
        Some(Command::Analyze(args)) => analyze::run(args),
        Some(Command::Completions { shell }) => {
            clap_complete::generate(shell, &mut Cli::command(), "rustafits", &mut std::io::stdout());
            Ok(true)
        }
    };
    match result {
        Ok(true) => {}
        Ok(false) => std::process::exit(1),
        Err(e) => {
            eprintln!("error: {e:#}");
            std::process::exit(1);
        }
    }
}
