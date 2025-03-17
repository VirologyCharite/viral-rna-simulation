import sys
import polars as pl
import plotly.express as px

from viral_rna_simulation.cells import Cells


def make_plot(cells: Cells, filename: str):
    positive_changes, negative_changes = cells.mutation_counts()
    overall_changes = positive_changes + negative_changes

    if overall_changes:
        # title = cells.summary().replace("\n", "<br>")

        from_positive, from_negative = cells.apparent_mutation_counts()
        apparent_changes = from_positive + from_negative
        assert apparent_changes

        TRANSITIONS = "AG", "GA", "CT", "TC"
        TRANSVERSIONS = "AT", "TA", "AC", "CA", "GT", "TG", "GC", "CG"

        BARCHART_CATEGORIES = [f"{t[0]}->{t[1]}" for t in TRANSITIONS] + [
            f"{t[0]}->{t[1]}" for t in TRANSVERSIONS
        ]
        changes = []
        origins = []
        counts = []

        for change in TRANSITIONS + TRANSVERSIONS:
            for origin, counter in (
                    ("Actual overall", overall_changes),
                    ("Actual (+) RNA", positive_changes),
                    ("Actual (-) RNA", negative_changes),
                    ("Apparent", apparent_changes),
            ):
                changes.append(f"{change[0]}->{change[1]}")
                origins.append(origin)
                counts.append(counter[change])

        df = pl.DataFrame({
            "Change": changes,
            "Origin": origins,
            "Count": counts,
        })

        fig = px.bar(
            df,
            x="Change",
            y="Count",
            color="Origin",
            barmode="group",
            category_orders={"Change": BARCHART_CATEGORIES},
            height=300,
            # title=title,
        )

        if filename.endswith(".html"):
            fig.write_html(filename)
        else:
            fig.write_image(filename)

        print(f"Wrote plot to {filename!r}.", file=sys.stderr)
    else:
        print(f"No genetic changes found, not writing {filename!r}.", file=sys.stderr)
