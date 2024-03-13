import run_tests
import os


commands = run_tests.generate_commands(
        "--algorithm concorde",
        os.path.join("test_results", "concorde"),
        "int(row['Number of vertices']) < 1000")


if __name__ == "__main__":
    for command in commands:
        run_tests.run(command)
