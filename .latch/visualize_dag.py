import json
import webbrowser
from pathlib import Path

import dagviz
import networkx as nx

with open(".latch/dag.json", "r") as f:
    dag = json.load(f)


G = nx.DiGraph()
for v in dag["vertices"]:
    c = v["content"]
    label = c["label"] if c["label"] is not None else c["type"]
    typ = c["type"]

    G.add_node(c["id"], label=f"{typ}: {label}")

for e in dag["edges"]:
    c = e["content"]
    src, dest = c["from"], c["to"]

    if src is None or dest is None:
        print(e)
        continue

    label = c["label"]
    attrs = {}
    if label is not None:
        attrs["label"] = label

    G.add_edge(src, dest, **attrs)

dag_path = Path(".latch/dag.svg").resolve()
with open(dag_path, "w") as f:
    f.write(dagviz.render_svg(G))

webbrowser.open(f"file://{dag_path}")
