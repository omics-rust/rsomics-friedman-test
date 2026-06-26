use std::fs::File;
use std::io::{self, BufReader};
use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use rsomics_common::{CommonFlags, RsomicsError, ToolMeta, run};

use rsomics_friedman_test::{friedman, parse_matrix};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

/// Friedman chi-square test for repeated measures (`scipy.stats.friedmanchisquare`).
///
/// Input is a TSV/whitespace matrix: one row per block (subject), one column per
/// treatment. At least three treatment columns are required. Values are ranked
/// within each block. Output is a single line `Q<TAB>df<TAB>p`.
#[derive(Parser, Debug)]
#[command(name = "rsomics-friedman-test", version, about, long_about = None)]
pub struct Cli {
    /// Input matrix (rows = blocks, cols = treatments); `-` or omitted reads stdin.
    #[arg(value_name = "DATA")]
    pub data: Option<PathBuf>,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn run(self) -> ExitCode {
        let common = self.common.clone();
        run(&common, META, || {
            let result = match &self.data {
                Some(p) if p.as_os_str() != "-" => {
                    let f = File::open(p).map_err(RsomicsError::Io)?;
                    let matrix = parse_matrix(BufReader::new(f))?;
                    friedman(&matrix)?
                }
                _ => {
                    let stdin = io::stdin();
                    let matrix = parse_matrix(stdin.lock())?;
                    friedman(&matrix)?
                }
            };
            if !common.json {
                println!("{}\t{}\t{}", result.q, result.df, result.p);
            }
            Ok(result)
        })
    }
}

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
