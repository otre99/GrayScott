use crate::gs_core::GrayScottSolver;
use clap::Parser;

mod gs_core;

/// Simple program to solve the
/// Gray-Scott problem  
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short, long, default_value_t = 0.055)]
    fparam: f64,

    #[arg(short, long, default_value_t = 0.062)]
    kparam: f64,

    #[arg(long, default_value_t = 0.16)]
    mu: f64,

    #[arg(long, default_value_t = 0.08)]
    mv: f64,

    #[arg(long, default_value_t = 256)]
    dim: usize,

    #[arg(long, default_value_t = 1000)]
    nsteps: usize,
}

fn main() {
    let args = Args::parse();

    let dim = args.dim;
    let mu = args.mu;
    let mv = args.mv;
    let fparam = args.fparam;
    let kparam = args.kparam;

    println!("Parameters: ");
    println!("  Dims   -> {}x{}", dim, dim);
    println!("  mu     -> {}", mu);
    println!("  mv     -> {}", mv);
    println!("  kparam -> {}", kparam);
    println!("  fparam -> {}", fparam);
    println!("  nspeps -> {}", args.nsteps);

    let mut gs = GrayScottSolver::new(dim, mu, mv, fparam, kparam);

    gs.save_solution("./");
    for _ in 0..args.nsteps {
        gs.step_seq();
    }
    gs.save_solution("./");
}
