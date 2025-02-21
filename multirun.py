import os
import sys
import argparse


def read_config(file):
    with open(file, "r") as f:
        config = f.read()

    config = config.split("\n")
    config = [line for line in config if not line.startswith("#")]

    config_dict = {}
    for line in config:
        line = [x.strip() for x in line.split("=")]

        if len(line) != 2:
            continue

        k = line[0]
        v = line[1]
        config_dict[k] = v

    return config_dict


def save_config(config, file):
    with open(file, "w") as f:
        for k, v in config.items():
            f.write(f"{k} = {v}\n")


def main():
    parser = argparse.ArgumentParser(description="Run flip at multiple pressures")
    parser.add_argument(
        "-p", "--pressures", type=float, nargs="+", help="Pressures to run flip at"
    )
    parser.add_argument(
        "-f",
        "--first",
        type=str,
        help="setup type for fist pressure (default: 'cubic)",
        default="cubic",
    )
    args = parser.parse_args()

    pressures = args.pressures
    first = args.first

    config = read_config("config.cfg")

    for i, p in enumerate(pressures):
        print(f"Running flip at pressure {p}")

        config["pressure"] = p
        if i != 0:
            prev = pressures[i - 1]
            config["setup"] = "restart"
            config["restartDir"] = f'"pressure-{str(prev).replace(".", "_")}"'
            print(f"Restarting from {config['restartDir']}")
        else:
            config["setup"] = f'"{first}"'

        config["dirName"] = f'"pressure-{str(p).replace(".", "_")}"'

        save_config(config, "config.cfg")

        os.system("./build/flip")


if __name__ == "__main__":
    main()
