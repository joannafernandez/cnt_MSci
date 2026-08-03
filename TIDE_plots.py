#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 13:35:51 2026

@author: jofernandez

Written with assistance from GPT 5.6 Thinking
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os 

#import data file and name cols 
def Find(file):
    genes = pd.read_csv(file, delimiter=",")
    genes.rename(columns={'Unnamed: 0': 'replicate'}, inplace=True)
    genes.rename(columns={'Unnamed: 1': 'indels'}, inplace=True)

    return genes

####s139a only
R0perc = Find('/Users/patricfernandez/ng_perc.csv')
R0perc= R0perc.dropna() #remove empty rows
R0pval = Find('/Users/patricfernandez/ng_pval.csv')
R0pval= R0pval.dropna() #remove empty rows

R0perc.loc[:, ~R0perc.columns.isin(['replicate', 'indels'])] = R0perc.loc[:, ~R0perc.columns.isin(['replicate', 'indels'])].mask(R0pval > 0.001, 0)


melted_df = R0perc.melt(id_vars=['replicate', 'indels'],
                               var_name='condition',
                               value_vars= [ 'GAPDH_T', 'GAPDH_puc', 'ACTIN_T',
                                      'ACTIN_puc'],
                               value_name='percentage')
melted_df['percentage'] = melted_df['percentage'].astype(int)

# Make sure indel sizes are numeric
melted_df["indels"] = pd.to_numeric(
    melted_df["indels"],
    errors="coerce"
)

# Assign repair pathway based on indel size
def classify_indel(indel):
    if pd.isna(indel):
        return None
    elif indel == 0:
        return "uncut"
    elif indel == 9:          # +9 indel
        return "HR"
    elif -4 <= indel <= 2:
        return "NHEJ"
    elif -20 <= indel <= -3:
        return "MMEJ"
    else:
        return None


melted_df["class"] = melted_df["indels"].apply(classify_indel)

# Remove indels that were not assigned to a repair class
melted_df = melted_df[
    melted_df["class"].notna()
].copy()

# Sum the percentages belonging to each repair class
summed_df = (
    melted_df
    .groupby(
        ["condition", "replicate", "class"],
        as_index=False
    )
    .agg(
        summed_percentage=("percentage", "sum")
    )
)

print(summed_df)

subset = summed_df.copy()

#Calculate total NHEJ + uncut per condition/replicate
total = subset.groupby(['condition', 'replicate'])['summed_percentage'].transform('sum')
#Calculate percentage within that total
subset['relative_percentage'] = (subset['summed_percentage'] / total) * 100
# Split 'condition' into two new columns: 'condition' and 'time'
subset[['condition', 'time']] = subset['condition'].str.split('_', n=1, expand=True)


#find total % of sequences that are either NHEJ or MMEJ. (uncut/wt sequences are only ever "0")
summed_df = (
    melted_df
    .groupby(['condition', 'replicate', 'class'], as_index=False)
    .agg(summed_percentage=('percentage', 'sum'))
)

subset = summed_df[summed_df['class'].isin(['MMEJ', 'NHEJ', 'uncut', "HR"])].copy()

#Calculate total NHEJ + uncut per condition/replicate
total = subset.groupby(['condition', 'replicate'])['summed_percentage'].transform('sum')
#Calculate percentage within that total
subset['relative_percentage'] = (subset['summed_percentage'] / total) * 100
# Split 'condition' into two new columns: 'condition' and 'time'
subset[['condition', 'time']] = subset['condition'].str.split('_', n=1, expand=True)


#%%

#----------------------------
# 1) Collapse to mutation vs uncut
# ----------------------------
df0 = subset.copy()


# ----------------------------
# 2) Categorical order + palette
# ----------------------------
time_order = ["T", "puc"]          # or your existing time_order
condition_order = ["ACTIN", "GAPDH"]        # or your existing condition_order
class_order = ['MMEJ', 'NHEJ', 'uncut', "HR"]    # NEW

df0["time"] = pd.Categorical(df0["time"], categories=time_order, ordered=True)
df0["condition"] = pd.Categorical(df0["condition"], categories=condition_order, ordered=True)
df0["class"] = pd.Categorical(df0["class"], categories=class_order, ordered=True)

# palette (use yours if you want)
custom_palette = {
    "MMEJ": "#DC6B83",
    "NHEJ": "#D4D3CF",
    "uncut": "#F3C773",
    "HR": "75B1CE"
}

# ----------------------------
# 3) Means for stacked bars
# ----------------------------
means = (
    df0.groupby(["time", "condition", "class"], observed=True)["relative_percentage"]
      .mean()
      .unstack("class")
      .reindex(index=pd.MultiIndex.from_product([time_order, condition_order],
                                                names=["time", "condition"]),
               columns=class_order)
      .fillna(0.0)
)


#%%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


# ============================================================
# 1. SUM INDELS INTO REPAIR CLASSES
# ============================================================

summed_df = (
    melted_df
    .groupby(
        ["condition", "replicate", "class"],
        as_index=False
    )
    .agg(
        summed_percentage=("percentage", "sum")
    )
)


# Keep the four classes used in the composition
subset = summed_df[
    summed_df["class"].isin(
        ["MMEJ", "NHEJ", "HR"]
    )
].copy()


# ============================================================
# 2. CALCULATE RELATIVE PERCENTAGES
# ============================================================

# Total selected percentage for each condition and replicate
subset["total_selected"] = (
    subset
    .groupby(
        ["condition", "replicate"]
    )["summed_percentage"]
    .transform("sum")
)

# Express each repair class as a percentage of that total
subset["relative_percentage"] = (
    subset["summed_percentage"]
    / subset["total_selected"]
    * 100
)




subset[["locus", "guide"]] = (
    subset["condition"]
    .str.split("_", n=1, expand=True)
)


# ============================================================
# 4. ORDERS AND COLOURS
# ============================================================

locus_order = ["ACTIN", "GAPDH"]

guide_order = {
    "ACTIN": ["T", "puc"],
    "GAPDH": ["T", "puc"]
}

# Stack order from bottom to top
class_order = [
    "NHEJ",
    "MMEJ",
    "HR",
    "uncut"
]

custom_palette = {
    "NHEJ": "#D4D3CF",
    "MMEJ": "#DC6B83",
    "HR": "#75B1CE",
    "uncut": "#F3C773"
}

rep_palette = {
    1: "#ed3658",
    2: "#36bbac",
    3: "#FFFFFF"
}


# Make replicate numeric
subset["replicate"] = pd.to_numeric(
    subset["replicate"],
    errors="raise"
).astype(int)




# ============================================================
# 5. ENSURE MISSING CLASSES ARE INCLUDED AS ZERO
# ============================================================

# This is important when one replicate has no indels belonging
# to one of the repair classes.

complete_rows = []

for locus in locus_order:
    for guide in guide_order[locus]:

        condition_name = f"{locus}_{guide}"

        condition_data = subset[
            subset["condition"] == condition_name
        ]

        replicates = sorted(
            condition_data["replicate"].unique()
        )

        for replicate in replicates:
            for repair_class in class_order:

                existing = condition_data[
                    (condition_data["replicate"] == replicate) &
                    (condition_data["class"] == repair_class)
                ]

                if existing.empty:
                    value = 0.0
                else:
                    value = existing[
                        "relative_percentage"
                    ].iloc[0]

                complete_rows.append({
                    "condition": condition_name,
                    "locus": locus,
                    "guide": guide,
                    "replicate": replicate,
                    "class": repair_class,
                    "relative_percentage": value
                })

df0 = pd.DataFrame(complete_rows)


# ============================================================
# 6. MEAN VALUES FOR STACKED BARS
# ============================================================

means = (
    df0
    .groupby(
        ["locus", "guide", "class"],
        as_index=False
    )["relative_percentage"]
    .mean()
)

#%% without uncut  

# ============================================================
# ORDERS AND COLOURS
# ============================================================

locus_order = ["ACTIN", "GAPDH"]

guide_order = {
    "ACTIN": ["T", "puc"],
    "GAPDH": ["T", "puc"]
}

class_order = [
    "NHEJ",
    "MMEJ",
    "HR"
]

custom_palette = {
    "NHEJ": "#D4D3CF",
    "MMEJ": "#DC6B83",
    "HR": "#75B1CE"
}

rep_palette = {
    1: "#ed3658",
    2: "#36bbac",
    3: "#FFFFFF"
}


# Make sure replicate is numeric
df0["replicate"] = pd.to_numeric(
    df0["replicate"],
    errors="raise"
).astype(int)


# ============================================================
# CALCULATE MEAN AND SEM
# ============================================================

summary = (
    df0
    .groupby(
        ["locus", "guide", "class"],
        as_index=False
    )
    .agg(
        mean_percentage=("relative_percentage", "mean"),
        sem_percentage=(
            "relative_percentage",
            lambda x: x.std(ddof=1) / np.sqrt(len(x))
            if len(x) > 1 else 0
        )
    )
)


# ============================================================
# PLOT
# ============================================================

fig, axes = plt.subplots(
    nrows=1,
    ncols=2,
    figsize=(9, 5),
    sharey=True,
    gridspec_kw={"wspace": 0.20}
)

bar_width = 0.18

# Position of each repair class around the guide centre
class_offsets = {
    "NHEJ": -1.5 * bar_width,
    "MMEJ": -0.5 * bar_width,
    "HR": 0.5 * bar_width,
    "uncut": 1.5 * bar_width
}

# Small offsets so replicate dots do not completely overlap
replicate_offsets = {
    1: -0.035,
    2: 0,
    3: 0.035
}


for ax, locus in zip(axes, locus_order):

    current_guides = guide_order[locus]
    guide_positions = np.arange(len(current_guides))

    for guide_x, guide in zip(
        guide_positions,
        current_guides
    ):

        for repair_class in class_order:

            bar_x = (
                guide_x
                + class_offsets[repair_class]
            )

            # ----------------------------------------------
            # Mean and SEM for this bar
            # ----------------------------------------------
            summary_row = summary[
                (summary["locus"] == locus) &
                (summary["guide"] == guide) &
                (summary["class"] == repair_class)
            ]

            if summary_row.empty:
                mean_value = 0
                sem_value = 0
            else:
                mean_value = float(
                    summary_row["mean_percentage"].iloc[0]
                )

                sem_value = float(
                    summary_row["sem_percentage"].iloc[0]
                )

            ax.bar(
                bar_x,
                mean_value,
                width=bar_width,
                color=custom_palette[repair_class],
                edgecolor="black",
                linewidth=1.4,
                alpha=0.85,
                zorder=2
            )

            ax.errorbar(
                bar_x,
                mean_value,
                yerr=sem_value,
                fmt="none",
                ecolor="black",
                elinewidth=1.5,
                capsize=4,
                capthick=1.5,
                zorder=4
            )

            # ----------------------------------------------
            # Replicate points
            # ----------------------------------------------
            replicate_data = df0[
                (df0["locus"] == locus) &
                (df0["guide"] == guide) &
                (df0["class"] == repair_class)
            ].copy()

            for _, row in replicate_data.iterrows():

                replicate = int(row["replicate"])

                dot_x = (
                    bar_x
                    + replicate_offsets.get(
                        replicate,
                        0
                    )
                )

                ax.scatter(
                    dot_x,
                    row["relative_percentage"],
                    s=90,
                    facecolor=rep_palette.get(
                        replicate,
                        "white"
                    ),
                    edgecolor="black",
                    linewidth=1.4,
                    zorder=5,
                    clip_on=False
                )

    # ========================================================
    # AXIS FORMATTING
    # ========================================================

    ax.set_xticks(guide_positions)

    ax.set_xticklabels(
        current_guides,
        fontsize=12
    )

    ax.set_title(
        locus,
        fontsize=15,
        pad=10
    )

    ax.set_xlim(
        -0.6,
        len(current_guides) - 0.4
    )

    ax.set_ylim(0, 105)

    ax.set_yticks(
        [0, 20, 40, 60, 80, 100]
    )

    ax.spines["left"].set_linewidth(2)
    ax.spines["bottom"].set_linewidth(2)

    ax.tick_params(
        axis="both",
        width=1.6,
        length=4,
        labelsize=12
    )

    ax.set_xlabel("")


axes[0].set_ylabel(
    "Relative percentage",
    fontsize=14
)


# ============================================================
# LEGENDS
# ============================================================

class_handles = [
    Patch(
        facecolor=custom_palette[repair_class],
        edgecolor="black",
        label=repair_class
    )
    for repair_class in class_order
]

replicate_handles = [
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="none",
        markerfacecolor=rep_palette[replicate],
        markeredgecolor="black",
        markeredgewidth=1.4,
        markersize=9,
        label=f"Replicate {replicate}"
    )
    for replicate in sorted(rep_palette)
]

class_legend = axes[1].legend(
    handles=class_handles,
    title="Repair class",
    frameon=False,
    fontsize=10,
    title_fontsize=11,
    loc="upper left",
    bbox_to_anchor=(1.03, 1.00)
)

axes[1].add_artist(class_legend)

axes[1].legend(
    handles=replicate_handles,
    title="Biological replicate",
    frameon=False,
    fontsize=10,
    title_fontsize=11,
    loc="upper left",
    bbox_to_anchor=(1.03, 0.55)
)

sns.despine(fig=fig)

plt.tight_layout()
plt.show()





#%%
#%%
#statistics




statsdf = subset.copy()


stats_wide = (
    statsdf
    .pivot_table(
        index=[
            "condition",
            "replicate",
            "locus",
            "guide"
        ],
        columns="class",
        values="relative_percentage",
        aggfunc="sum",
        fill_value=0
    )
    .reset_index()
)

# Remove the extra name attached to the columns index
stats_wide.columns.name = None

print(stats_wide)

#statsdf['lir'] = np.log2(statsdf['NHEJ']/(statsdf['MMEJ']*statsdf['HR']))


pseudocount = 0.01

stats_wide["NHEJ_balance"] = np.log2(
    (stats_wide["NHEJ"] + pseudocount)
    / np.sqrt(
        (stats_wide["MMEJ"] + pseudocount)
        * (stats_wide["HR"] + pseudocount)
    )
)


#%%


from scipy import stats
import pandas as pd
import numpy as np

# Use whichever column name you created:
ilr_col = "NHEJ_balance"
# Or, for your unscaled log2 balance:
# ilr_col = "NHEJ_balance"


# Keep only the two targeted DSB conditions
test_df = stats_wide[
    stats_wide["condition"].isin([
        "ACTIN_T",
        "GAPDH_T"
    ])
].copy()


# Put the two loci into separate columns, matched by replicate
paired_df = (
    test_df
    .pivot(
        index="replicate",
        columns="condition",
        values=ilr_col
    )
    .dropna(
        subset=["ACTIN_T", "GAPDH_T"]
    )
)

print(paired_df)

#%%


ttest = stats.ttest_rel(
    paired_df["GAPDH_T"],
    paired_df["ACTIN_T"]
)

print(f"Paired t-test: t = {ttest.statistic:.3f}")
print(f"p = {ttest.pvalue:.4f}")


#%%

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D

# ============================================================
# COLOURS
# ============================================================

condition_palette = {
    "ACTIN_T": "#B3B3B3",
    "GAPDH_T": "#75B1CE"
}

rep_palette = {
    1: "#ed3658",
    2: "#36bbac",
    3: "#FFFFFF"
}

condition_order = [
    "ACTIN_T",
    "GAPDH_T"
]

x_positions = {
    "ACTIN_T": 0,
    "GAPDH_T": 1
}


# ============================================================
# PLOT
# ============================================================

fig, ax = plt.subplots(
    figsize=(4.5, 5.2)
)

# Plot paired lines first
for replicate, row in paired_df.iterrows():

    ax.plot(
        [x_positions["ACTIN_T"], x_positions["GAPDH_T"]],
        [row["ACTIN_T"], row["GAPDH_T"]],
        color="grey",
        linewidth=1.4,
        alpha=0.7,
        zorder=1
    )

    # Replicate dots
    for condition in condition_order:

        ax.scatter(
            x_positions[condition],
            row[condition],
            s=150,
            facecolor=rep_palette[int(replicate)],
            edgecolor="black",
            linewidth=1.6,
            zorder=4
        )


# ============================================================
# MEAN ± SEM
# ============================================================

for condition in condition_order:

    values = paired_df[condition].to_numpy(dtype=float)

    mean_value = values.mean()
    sem_value = (
        values.std(ddof=1)
        / np.sqrt(len(values))
    )

    # Mean horizontal line
    ax.plot(
        [
            x_positions[condition] - 0.17,
            x_positions[condition] + 0.17
        ],
        [mean_value, mean_value],
        color=condition_palette[condition],
        linewidth=4,
        zorder=5
    )

    # SEM
    ax.errorbar(
        x_positions[condition],
        mean_value,
        yerr=sem_value,
        fmt="none",
        ecolor=condition_palette[condition],
        elinewidth=2,
        capsize=6,
        capthick=2,
        zorder=5
    )


# ============================================================
# FORMATTING
# ============================================================

ax.set_xticks([0, 1])

ax.set_xticklabels(
    ["ACTIN", "GAPDH"],
    fontsize=13
)

ax.set_ylabel(
    r"$\log_{2}$ NHEJ / geometric mean(MMEJ, HR)",
    fontsize=13
)

ax.set_xlim(-0.45, 1.45)

# Set automatically from data, with some extra space
data_min = paired_df.min().min()
data_max = paired_df.max().max()
padding = (data_max - data_min) * 0.15

ax.set_ylim(
    data_min - padding,
    data_max + padding
)

ax.spines["left"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)

ax.tick_params(
    axis="both",
    width=1.6,
    length=4,
    labelsize=12
)


# ============================================================
# REPLICATE LEGEND
# ============================================================

replicate_handles = [
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="none",
        markerfacecolor=rep_palette[replicate],
        markeredgecolor="black",
        markeredgewidth=1.4,
        markersize=9,
        label=f"Replicate {replicate}"
    )
    for replicate in sorted(rep_palette)
]

ax.legend(
    handles=replicate_handles,
    title="Biological replicate",
    frameon=False,
    fontsize=10,
    title_fontsize=11,
    loc="upper left",
    bbox_to_anchor=(1.02, 1)
)

sns.despine()

plt.tight_layout()
plt.show()





#%%
#for plotting only, with uncut

# ============================================================
# 1. SUM INDELS INTO REPAIR CLASSES
# ============================================================

summed_df = (
    melted_df
    .groupby(
        ["condition", "replicate", "class"],
        as_index=False
    )
    .agg(
        summed_percentage=("percentage", "sum")
    )
)


# Keep the four classes used in the composition
subset = summed_df[
    summed_df["class"].isin(
        ["MMEJ", "NHEJ", "HR", "uncut"]
    )
].copy()


# ============================================================
# 2. CALCULATE RELATIVE PERCENTAGES
# ============================================================

# Total selected percentage for each condition and replicate
subset["total_selected"] = (
    subset
    .groupby(
        ["condition", "replicate"]
    )["summed_percentage"]
    .transform("sum")
)

# Express each repair class as a percentage of that total
subset["relative_percentage"] = (
    subset["summed_percentage"]
    / subset["total_selected"]
    * 100
)




subset[["locus", "guide"]] = (
    subset["condition"]
    .str.split("_", n=1, expand=True)
)


# ============================================================
# 4. ORDERS AND COLOURS
# ============================================================

locus_order = ["ACTIN", "GAPDH"]

guide_order = {
    "ACTIN": ["T", "puc"],
    "GAPDH": ["T", "puc"]
}

# Stack order from bottom to top
class_order = [
    "NHEJ",
    "MMEJ",
    "HR",
    "uncut"
]

custom_palette = {
    "NHEJ": "#D4D3CF",
    "MMEJ": "#DC6B83",
    "HR": "#75B1CE",
    "uncut": "#F3C773"
}

rep_palette = {
    1: "#ed3658",
    2: "#36bbac",
    3: "#FFFFFF"
}


# Make replicate numeric
subset["replicate"] = pd.to_numeric(
    subset["replicate"],
    errors="raise"
).astype(int)




# ============================================================
# 5. ENSURE MISSING CLASSES ARE INCLUDED AS ZERO
# ============================================================

# This is important when one replicate has no indels belonging
# to one of the repair classes.

complete_rows = []

for locus in locus_order:
    for guide in guide_order[locus]:

        condition_name = f"{locus}_{guide}"

        condition_data = subset[
            subset["condition"] == condition_name
        ]

        replicates = sorted(
            condition_data["replicate"].unique()
        )

        for replicate in replicates:
            for repair_class in class_order:

                existing = condition_data[
                    (condition_data["replicate"] == replicate) &
                    (condition_data["class"] == repair_class)
                ]

                if existing.empty:
                    value = 0.0
                else:
                    value = existing[
                        "relative_percentage"
                    ].iloc[0]

                complete_rows.append({
                    "condition": condition_name,
                    "locus": locus,
                    "guide": guide,
                    "replicate": replicate,
                    "class": repair_class,
                    "relative_percentage": value
                })

df0 = pd.DataFrame(complete_rows)


# ============================================================
# 6. MEAN VALUES FOR STACKED BARS
# ============================================================

means = (
    df0
    .groupby(
        ["locus", "guide", "class"],
        as_index=False
    )["relative_percentage"]
    .mean()
)

# ============================================================
# ORDERS AND COLOURS
# ============================================================

locus_order = ["ACTIN", "GAPDH"]

guide_order = {
    "ACTIN": ["T", "puc"],
    "GAPDH": ["T", "puc"]
}

class_order = [
    "NHEJ",
    "MMEJ",
    "HR",
    "uncut"
]

custom_palette = {
    "NHEJ": "#D4D3CF",
    "MMEJ": "#DC6B83",
    "HR": "#75B1CE",
    "uncut": "#F3C773"
}

rep_palette = {
    1: "#ed3658",
    2: "#36bbac",
    3: "#FFFFFF"
}


# Make sure replicate is numeric
df0["replicate"] = pd.to_numeric(
    df0["replicate"],
    errors="raise"
).astype(int)


# ============================================================
# CALCULATE MEAN AND SEM
# ============================================================

summary = (
    df0
    .groupby(
        ["locus", "guide", "class"],
        as_index=False
    )
    .agg(
        mean_percentage=("relative_percentage", "mean"),
        sem_percentage=(
            "relative_percentage",
            lambda x: x.std(ddof=1) / np.sqrt(len(x))
            if len(x) > 1 else 0
        )
    )
)


# ============================================================
# PLOT
# ============================================================

fig, axes = plt.subplots(
    nrows=1,
    ncols=2,
    figsize=(9, 5),
    sharey=True,
    gridspec_kw={"wspace": 0.20}
)

bar_width = 0.18

# Position of each repair class around the guide centre
class_offsets = {
    "NHEJ": -1.5 * bar_width,
    "MMEJ": -0.5 * bar_width,
    "HR": 0.5 * bar_width,
    "uncut": 1.5 * bar_width
}

# Small offsets so replicate dots do not completely overlap
replicate_offsets = {
    1: -0.035,
    2: 0,
    3: 0.035
}


for ax, locus in zip(axes, locus_order):

    current_guides = guide_order[locus]
    guide_positions = np.arange(len(current_guides))

    for guide_x, guide in zip(
        guide_positions,
        current_guides
    ):

        for repair_class in class_order:

            bar_x = (
                guide_x
                + class_offsets[repair_class]
            )

            # ----------------------------------------------
            # Mean and SEM for this bar
            # ----------------------------------------------
            summary_row = summary[
                (summary["locus"] == locus) &
                (summary["guide"] == guide) &
                (summary["class"] == repair_class)
            ]

            if summary_row.empty:
                mean_value = 0
                sem_value = 0
            else:
                mean_value = float(
                    summary_row["mean_percentage"].iloc[0]
                )

                sem_value = float(
                    summary_row["sem_percentage"].iloc[0]
                )

            ax.bar(
                bar_x,
                mean_value,
                width=bar_width,
                color=custom_palette[repair_class],
                edgecolor="black",
                linewidth=1.4,
                alpha=0.85,
                zorder=2
            )

            ax.errorbar(
                bar_x,
                mean_value,
                yerr=sem_value,
                fmt="none",
                ecolor="black",
                elinewidth=1.5,
                capsize=4,
                capthick=1.5,
                zorder=4
            )

            # ----------------------------------------------
            # Replicate points
            # ----------------------------------------------
            replicate_data = df0[
                (df0["locus"] == locus) &
                (df0["guide"] == guide) &
                (df0["class"] == repair_class)
            ].copy()

            for _, row in replicate_data.iterrows():

                replicate = int(row["replicate"])

                dot_x = (
                    bar_x
                    + replicate_offsets.get(
                        replicate,
                        0
                    )
                )

                ax.scatter(
                    dot_x,
                    row["relative_percentage"],
                    s=90,
                    facecolor=rep_palette.get(
                        replicate,
                        "white"
                    ),
                    edgecolor="black",
                    linewidth=1.4,
                    zorder=5,
                    clip_on=False
                )

    # ========================================================
    # AXIS FORMATTING
    # ========================================================

    ax.set_xticks(guide_positions)

    ax.set_xticklabels(
        current_guides,
        fontsize=12
    )

    ax.set_title(
        locus,
        fontsize=15,
        pad=10
    )

    ax.set_xlim(
        -0.6,
        len(current_guides) - 0.4
    )

    ax.set_ylim(0, 105)

    ax.set_yticks(
        [0, 20, 40, 60, 80, 100]
    )

    ax.spines["left"].set_linewidth(2)
    ax.spines["bottom"].set_linewidth(2)

    ax.tick_params(
        axis="both",
        width=1.6,
        length=4,
        labelsize=12
    )

    ax.set_xlabel("")


axes[0].set_ylabel(
    "Relative percentage",
    fontsize=14
)


# ============================================================
# LEGENDS
# ============================================================

class_handles = [
    Patch(
        facecolor=custom_palette[repair_class],
        edgecolor="black",
        label=repair_class
    )
    for repair_class in class_order
]

replicate_handles = [
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="none",
        markerfacecolor=rep_palette[replicate],
        markeredgecolor="black",
        markeredgewidth=1.4,
        markersize=9,
        label=f"Replicate {replicate}"
    )
    for replicate in sorted(rep_palette)
]

class_legend = axes[1].legend(
    handles=class_handles,
    title="Repair class",
    frameon=False,
    fontsize=10,
    title_fontsize=11,
    loc="upper left",
    bbox_to_anchor=(1.03, 1.00)
)

axes[1].add_artist(class_legend)

axes[1].legend(
    handles=replicate_handles,
    title="Biological replicate",
    frameon=False,
    fontsize=10,
    title_fontsize=11,
    loc="upper left",
    bbox_to_anchor=(1.03, 0.55)
)

sns.despine(fig=fig)

plt.tight_layout()
plt.show()
