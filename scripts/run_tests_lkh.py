import run_tests
import os


commands = run_tests.generate_commands(
        "--algorithm lkh"
        " --candidate-set-type POPMUSIC"
        " --initial-period 100"
        " --runs 1"
        " --max-trials 1",
        os.path.join("test_results", "lkh"))


if __name__ == "__main__":
    for command in commands:
        run_tests.run(command)
