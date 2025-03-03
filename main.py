import numpy as np
from tqdm.auto import tqdm
from itertools import combinations, product
from collections import OrderedDict
import pandas as pd
import multiprocessing as mp
import os
import fire
import pyrootutils

PROJECT_ROOT = path = pyrootutils.find_root(
    search_from=__file__, indicator=".project-root"
)

NUM_ELEMENTS_IN_ASSAY = 7
NUM_VALUES_PER_ELEMENT = 20

INTERACTION_TABLE_SIZE = NUM_VALUES_PER_ELEMENT * NUM_ELEMENTS_IN_ASSAY
NUM_INTERACTIONS = NUM_VALUES_PER_ELEMENT**NUM_ELEMENTS_IN_ASSAY

INPUT_DATA_FILENAME = "clean_data.csv"
PROCESSED_INTERACTIONS_FILENAME = "scores.npy"
TARGET_QUALITY_FILENAME = "target_quality.tsv"


def to_int(elements: list) -> int:
    assert len(elements) == NUM_ELEMENTS_IN_ASSAY
    assert all(0 <= x < NUM_VALUES_PER_ELEMENT for x in elements)
    return sum(x * NUM_VALUES_PER_ELEMENT**i for i, x in enumerate(elements))


def to_list(value: int) -> list:
    assert 0 <= value < NUM_INTERACTIONS
    return [
        (value // NUM_VALUES_PER_ELEMENT**i) % NUM_VALUES_PER_ELEMENT
        for i in range(NUM_ELEMENTS_IN_ASSAY)
    ]


def to_scaled_indices(i: int) -> list:
    elements = to_list(i)
    return [val + i * NUM_VALUES_PER_ELEMENT for i, val in enumerate(elements)]


# def scaled_indices_to_int(scaled_indices: list) -> int:



def get_score(i: int, interaction_table: np.ndarray) -> tuple:
    scaled_indices = to_scaled_indices(i)
    interaction_indices = list(combinations(scaled_indices, 2))    
    rows, cols = zip(*interaction_indices)
    scores = interaction_table[rows, cols]
    total = np.sum(scores)
    return (i, total)


def process(parallel: bool = False):
    interaction_table_raw = pd.read_csv(INPUT_DATA_FILENAME)
    
    interaction_table = (
        interaction_table_raw.drop(columns=["Unnamed: 0", "SUM"])
        .apply(pd.to_numeric)
        .to_numpy()
    )

    interaction_table[np.isnan(interaction_table)] = 0
    interaction_table_sym = (
        interaction_table + interaction_table.T - np.diag(interaction_table.diagonal())
    )

    assert interaction_table_sym.shape == (INTERACTION_TABLE_SIZE, INTERACTION_TABLE_SIZE)

    get_score_fn = lambda i: get_score(i, interaction_table_sym)

    if parallel:
        num_processes = max(2, os.cpu_count() - 2)
        with mp.Pool(num_processes) as pool:
            scores = np.array(
                list(
                    tqdm(
                        pool.imap(get_score_fn, range(NUM_INTERACTIONS)),
                        total=NUM_INTERACTIONS,
                    )
                )
            )
        with open(PROCESSED_INTERACTIONS_FILENAME, "wb") as f:
            np.save(f, scores)
    else:
        scores = np.zeros(NUM_INTERACTIONS, dtype=np.int16)
        for i in tqdm(range(NUM_INTERACTIONS)):
            _, scores[i] = get_score_fn(i)

        with open(PROJECT_ROOT / PROCESSED_INTERACTIONS_FILENAME, "wb") as f:
            np.save(f, scores)


def analyze(
    scores_file: str = PROCESSED_INTERACTIONS_FILENAME,
    target_file: str = TARGET_QUALITY_FILENAME,
):

    # Load target quality data
    target_quality_df = pd.read_csv(PROJECT_ROOT / target_file, sep="\t")
    target_quality_df['good'] = target_quality_df['good'].astype(str)
    target_quality_df['bad'] = target_quality_df['bad'].astype(str)
    target_quality_df['no'] = target_quality_df['no'].astype(str)

    target_quality_df["good_proc"] = target_quality_df.apply(lambda x: [f'{i} ' + x['target'] for i in x['good'].split(', ')], axis=1)
    target_quality_df["bad_proc"] = target_quality_df.apply(lambda x: [f'{i} ' + x['target'] for i in x['bad'].split(', ')], axis=1)
    target_quality_df["no_proc"] = target_quality_df.apply(lambda x: [f'{i} ' + x['target'] for i in x['no'].split(', ')], axis=1)

    target_quality_df = (
        target_quality_df
            .drop(columns=["good", "bad", "no", "target"])
            .rename(columns={"good_proc": "good", "bad_proc": "bad", "no_proc": "no"})
            .melt(
                value_vars=["good", "bad", "no"], 
                var_name="quality", 
                value_name="design"
            )
            .explode("design")
            .reset_index(drop=True)[["design", "quality"]]
    )
    target_quality_df = target_quality_df[
        ~target_quality_df["design"].str.contains("nan ")
    ].reset_index(drop=True)

    # Filter for good designs
    good_designs_df = target_quality_df[target_quality_df["quality"] == "good"]
    good_designs_df["design"] = good_designs_df["design"].str.strip()
    oligose_table = pd.read_csv(INPUT_DATA_FILENAME)[["Unnamed: 0"]].rename(
        columns={"Unnamed: 0": "Oligose"}
    )["Oligose"].to_dict()
    oligose_table = {v: k for k, v in oligose_table.items()}
    good_designs_df["index"] = good_designs_df["design"].map(lambda x: oligose_table[x])
    good_simplex_indices = good_designs_df["index"].to_list()

    good_designs_df["target"] = good_designs_df["design"].str.split(" ").str[1]
    good_designs_df["design"] = good_designs_df["design"].str.split(" ").str[0].apply(pd.to_numeric)

    good_designs_dict = OrderedDict()
    for i, tgt in enumerate(["Av", "BVAB", "Catdp", "Cg", "Ck", "Meg", "Tv"]):
        good_designs_dict[i] = (tgt, good_designs_df[good_designs_df["target"] == tgt]["design"].to_list())
    good_designs_idx_lists = [good_designs_dict[i][1] for i in good_designs_dict.keys()]
    good_designs_idx_lists = [list(map(lambda x: x - 1, lst)) for lst in good_designs_idx_lists]
    good_interaction_combs = list(product(*good_designs_idx_lists))
    good_interaction_indices = [to_int(comb) for comb in good_interaction_combs]

    scores = np.load(PROJECT_ROOT / scores_file)
    print(f"Max: {np.max(scores)}")
    print(f"Zeroes: {np.sum(scores == 0)} ({np.sum(scores == 0) / len(scores):.2%})")

    top_scores_idx = np.argpartition(scores, -5)[-5:]
    top_scores = scores[top_scores_idx]
    print(top_scores, top_scores_idx)

    top_scores_zipped = list(zip(top_scores_idx, top_scores))
    top_scores_df = pd.DataFrame(top_scores_zipped, columns=["index", "score"])
    top_scores_df["list"] = top_scores_df["index"].apply(to_list).apply(lambda x: "; ".join(map(str, x)))
    print(top_scores_df)

    # top_scores_df.to_csv("highest_scores.csv", index=False, sep=",")
    top_scores_df.to_excel("highest_scores.xlsx", index=False)

    # all_interaction_scores_df = pd.DataFrame(scores, columns=["score"]).sort_values(by="score", ascending=False).head()
    # # all_interaction_scores_df["index"] = all_interaction_scores_df.index
    # # all_interaction_scores_df["list"] = all_interaction_scores_df["index"].apply(to_list)
    # print(all_interaction_scores_df)

    # raise SystemExit

    good_interaction_scores = []
    for idx in good_interaction_indices:
        good_interaction_scores.append({"index": idx, "score": scores[idx]})

    gis_df = pd.DataFrame(good_interaction_scores)
    gis_df["list"] = gis_df["index"].apply(to_list).apply(lambda x: "; ".join(map(str, x)))
    gis_df = gis_df.sort_values(by="score", ascending=True)
    print(gis_df.head())
    print(gis_df.tail())
    print(gis_df.shape)

    # gis_df.head().to_csv("allgood_lowest_scores.csv", index=False, sep=",")
    # gis_df.tail().to_csv("allgood_highest_scores.csv", index=False, sep=",")
    gis_df.head().to_excel("allgood_lowest_scores.xlsx", index=False)
    gis_df.tail().to_excel("allgood_highest_scores.xlsx", index=False)

    gis_zeros_df = gis_df[gis_df["score"] == 0]
    print(gis_zeros_df.shape)

    # gis_zeros_df.to_csv("allgood_zero_scores.csv", index=False, sep="\t")
    gis_zeros_df.to_excel("allgood_zero_scores.xlsx", index=False)

    # # sort good_interaction_scores by value
    # good_interaction_scores = dict(sorted(good_interaction_scores.items(), key=lambda x: x[1], reverse=True))

    # print(good_interaction_scores)

if __name__ == "__main__":
    fire.Fire()
