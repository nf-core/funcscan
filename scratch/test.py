import glob
import json
import os
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List, NamedTuple, Optional

from latch_cli.extras.nextflow.channel import get_boolean_value, get_mapper_outputs

# channel_dir = Path(".latch_compiled_channels")
# if channel_dir.exists():
#     shutil.rmtree(channel_dir)

# # os.environ["LATCH_TARGET_PROCESS_NAME"] = "testBash"
# os.environ["LATCH_PARAM_VALS"] = json.dumps(
#     [
#         [{"string": f"{chr(i)}1"} for i in range(65, 123)],
#         # [{"string": f"{chr(i)}2"} for i in range(65, 123)],
#     ]
# )
# # # os.environ["LATCH_PARAM_VALS"] = json.dumps([{"string": f"{chr(i)}"} for i in range(65, 123)])

# # os.environ["LATCH_EXPRESSION"] = (
# #     '{"ExpressionStatement":{"expression":{"BinaryExpression":{"leftExpression":{"VariableExpression":"res"},"operation":"=","rightExpression":{"MethodCallExpression":{"objectExpression":{"MethodCallExpression":{"objectExpression":{"VariableExpression":"Channel"},"method":"of","arguments":{"ArgumentListExpression":{"expressions":[]}}}},"method":"multiMap","arguments":{"ArgumentListExpression":{"expressions":[{"ClosureExpression":{"code":{"BlockStatement":{"statements":[{"ExpressionStatement":{"expression":{"MethodCallExpression":{"objectExpression":{"VariableExpression":"it"},"method":"toUpperCase","arguments":{"ArgumentListExpression":{"expressions":[]}}}},"labels":["upper"]}},{"ExpressionStatement":{"expression":{"MethodCallExpression":{"objectExpression":{"VariableExpression":"it"},"method":"toLowerCase","arguments":{"ArgumentListExpression":{"expressions":[]}}}},"labels":["lower"]}}],"scope":{"declaredVariables":[],"referencedClassVariables":[]},"labels":[]}},"parameters":[]}}]}}}}}},"labels":[]}}'
# # )
# # os.environ["LATCH_RETURN"] = r"""[
# #     "{\"ExpressionStatement\":{\"expression\":{\"BinaryExpression\":{\"leftExpression\":{\"VariableExpression\":\"upper\"},\"operation\":\"=\",\"rightExpression\":{\"PropertyExpression\":{\"objectExpression\":{\"VariableExpression\":\"res\"},\"property\":\"upper\"}}}},\"labels\":[]}}",
# #     "{\"ExpressionStatement\":{\"expression\":{\"BinaryExpression\":{\"leftExpression\":{\"VariableExpression\":\"lower\"},\"operation\":\"=\",\"rightExpression\":{\"PropertyExpression\":{\"objectExpression\":{\"VariableExpression\":\"res\"},\"property\":\"lower\"}}}},\"labels\":[]}}"
# # ]"""


pkg_root = Path.cwd()
# nf_script = pkg_root / "subworkflows" / "local" / "arg.nf"
# nf_script = pkg_root / "foo.nf"
nf_script = pkg_root / "main.nf"
# nf_script = pkg_root / "workflows" / "funcscan.nf"

# subprocess.run(
#     shlex.join(
#         [
#             str(pkg_root / ".latch/bin/nextflow"),
#             "run",
#             str(nf_script),
#             "-latchRegister",
#             # "-profile",
#             # "mamba",
#             # "-dump-channels",
#             # "--input",
#             # str(pkg_root / "assets" / "samplesheet.csv"),
#             # "--outdir",
#             # str(pkg_root),
#         ]
#     ),
#     executable="zsh",
#     shell=True,
# )


# # =========


channel_vals = [[
    {
        "list": [
            {"map": [{"key": {"string": "id"}, "value": {"string": "sample_1"}}]},
            {"path": "/root/work/1e/cc185813c2aa081041b844faacac31/wastewater_metagenome_contigs_1.fasta"},
            {"string": "23263"},
        ]
    },
    {
        "list": [
            {"map": [{"key": {"string": "id"}, "value": {"string": "sample_2"}}]},
            {"path": "/root/work/ae/1a4453e684560489f1f5a56e4c1e73/wastewater_metagenome_contigs_2.fasta"},
            {"string": "23263"},
        ]
    },
]]
# channel_vals = json.loads('[[{"string": "a"}, {"string": "b"}, {"string": "C"}]]')
# channel_vals = [[{"integer": 2}, {"integer": 3}]]

# channel_vals = [[{"path": "/root/assets/multiqc_config.yml"}]]
# channel_vals = [[{"string": "/root/assets/multiqc_config.yml"}]]
# channel_vals = [{"value": {"boolean": False}}, {"value": {"boolean": False}}]

channel_vals = [
    {"value": {"list": [{"path": "work/software_versions_mqc.yml"}]}},
    {"value": {"list": [{"path": "assets/multiqc_config.yml"}]}},
    {"value": {"list": []}},
    {"value": {"list": [{"path": "docs/images/nf-core-funcscan_logo_flat_light.png"}]}},
]

subprocess.run(
    [
        ".latch/bin/nextflow",
        "run",
        str(nf_script),
        "--input",
        str(pkg_root / "assets" / "samplesheet.csv"),
        "--outdir",
        str(pkg_root),
        "-profile",
        "mamba",
    ],
    env={
        **os.environ,
        "LATCH_INCLUDE_META": '{"path": "./modules/nf-core/multiqc/main.nf", "alias": "MULTIQC", "name": "MULTIQC"}',
        "LATCH_EXPRESSION": '{"ExpressionStatement":{"expression":{"MethodCallExpression":{"objectExpression":{"VariableExpression":"this"},"method":"MULTIQC","arguments":{"ArgumentListExpression":{"expressions":[{"MethodCallExpression":{"objectExpression":{"VariableExpression":"Channel"},"method":"placeholder","arguments":{"ArgumentListExpression":{"expressions":[]}}}},{"MethodCallExpression":{"objectExpression":{"VariableExpression":"Channel"},"method":"placeholder","arguments":{"ArgumentListExpression":{"expressions":[]}}}},{"MethodCallExpression":{"objectExpression":{"VariableExpression":"Channel"},"method":"placeholder","arguments":{"ArgumentListExpression":{"expressions":[]}}}},{"MethodCallExpression":{"objectExpression":{"VariableExpression":"Channel"},"method":"placeholder","arguments":{"ArgumentListExpression":{"expressions":[]}}}}]}}}},"labels":[]}}',
        "LATCH_RETURN": (
            '["{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"report\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"MULTIQC\\"},\\"property\\":\\"out\\"}},\\"operation\\":\\"[\\",\\"rightExpression\\":{\\"ConstantExpression\\":0}}}}},\\"labels\\":[]}}",'
            ' "{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"data\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"MULTIQC\\"},\\"property\\":\\"out\\"}},\\"operation\\":\\"[\\",\\"rightExpression\\":{\\"ConstantExpression\\":1}}}}},\\"labels\\":[]}}",'
            ' "{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"plots\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"MULTIQC\\"},\\"property\\":\\"out\\"}},\\"operation\\":\\"[\\",\\"rightExpression\\":{\\"ConstantExpression\\":2}}}}},\\"labels\\":[]}}",'
            ' "{\\"ExpressionStatement\\":{\\"expression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"VariableExpression\\":\\"versions\\"},\\"operation\\":\\"=\\",\\"rightExpression\\":{\\"BinaryExpression\\":{\\"leftExpression\\":{\\"PropertyExpression\\":{\\"objectExpression\\":{\\"VariableExpression\\":\\"MULTIQC\\"},\\"property\\":\\"out\\"}},\\"operation\\":\\"[\\",\\"rightExpression\\":{\\"ConstantExpression\\":3}}}}},\\"labels\\":[]}}"]'
        ),
        "LATCH_PARAM_VALS": json.dumps(channel_vals),
    },
    check=True,
)

out_channels = {}
files = [Path(f) for f in glob.glob(".latch/task-outputs/*.json")]

for file in files:
    out_channels[file.stem] = file.read_text()


shutil.rmtree(".latch/task-outputs")


class Respost_adapter_MULTIQC_1493_post(NamedTuple):
    report: Optional[str]
    data: Optional[str]
    plots: Optional[str]
    versions: Optional[str]


@dataclass
class Dataclass_1493_post:
    report: str
    data: str
    plots: str
    versions: str


default = [
    Dataclass_1493_post(
        report=out_channels.get(f"report", ""),
        data=out_channels.get(f"data", ""),
        plots=out_channels.get(f"plots", ""),
        versions=out_channels.get(f"versions", ""),
    )
]


output = get_mapper_outputs(Respost_adapter_MULTIQC_1493_post, default)
print(output)
