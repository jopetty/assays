import fire
import pyrootutils
import numpy as np


PROJECT_ROOT = path = pyrootutils.find_root(
    search_from=__file__, indicator=".project-root"
)


def main(
    data_file: str = "scores.npy"
):
    data_path = PROJECT_ROOT / data_file

    scores = np.load(data_path)
    print(scores[:10])

    num_values = len(scores)
    print(f"Number of values: {num_values}")

    max_value = np.max(scores)
    print(f"Max value: {max_value}")

    num_zeros = np.sum(scores == 0)
    print(f"Number of zeros: {num_zeros} ({num_zeros / num_values:.2%})")


if __name__ == "__main__":
    fire.Fire(main)