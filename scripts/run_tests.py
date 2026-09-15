import argparse
import os
import sys


parser = argparse.ArgumentParser(description='')
parser.add_argument('directory')
parser.add_argument(
        "-t", "--tests",
        type=str,
        nargs='*',
        help='')

args = parser.parse_args()


travelingsalesmansolver_main = os.path.join(
        "install",
        "bin",
        "travelingsalesmansolver")
data_dir = os.path.join("data", "tsplib")


# Small instances, so that both an exact solver (Concorde) and a
# heuristic (LKH) solve them quickly.
instances = [
        "ulysses16",
        "ulysses22",
        "att48",
        "berlin52",
        "eil51",
        "pr76",
        "rd100",
        "lin105",
        "kroC100",
        ]


if args.tests is None or "concorde" in args.tests:
    print("Concorde")
    print("--------")
    print()

    for instance in instances:
        instance_path = os.path.join(data_dir, instance + ".tsp")

        json_output_path = os.path.join(
                args.directory,
                "concorde",
                instance + ".json")
        if not os.path.exists(os.path.dirname(json_output_path)):
            os.makedirs(os.path.dirname(json_output_path))
        command = (
                travelingsalesmansolver_main
                + "  --verbosity-level 1"
                + "  --input \"" + instance_path + "\""
                + "  --algorithm \"concorde\""
                + "  --output \"" + json_output_path + "\"")
        print(command)
        status = os.system(command)
        if status != 0:
            sys.exit(1)
        print()
    print()
    print()


if args.tests is None or "lkh" in args.tests:
    print("LKH")
    print("---")
    print()

    for instance in instances:
        instance_path = os.path.join(data_dir, instance + ".tsp")

        json_output_path = os.path.join(
                args.directory,
                "lkh",
                instance + ".json")
        if not os.path.exists(os.path.dirname(json_output_path)):
            os.makedirs(os.path.dirname(json_output_path))
        command = (
                travelingsalesmansolver_main
                + "  --verbosity-level 1"
                + "  --input \"" + instance_path + "\""
                + "  --algorithm \"lkh\""
                + "  --time-limit 10"
                + "  --runs 1"
                + "  --output \"" + json_output_path + "\"")
        print(command)
        status = os.system(command)
        if status != 0:
            sys.exit(1)
        print()
    print()
    print()
