from pathlib import Path
from typing import Optional

from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import (
    LatchAuthor,
    NextflowFileParameter,
    NextflowMetadata,
    NextflowParameter,
)

NextflowMetadata(
    display_name="nf-core/funcscan",
    author=LatchAuthor(),
    parameters={
        "input": NextflowFileParameter(
            type=LatchFile,
            default=LatchFile("latch://1721.account/funcscan.samplesheet.csv"),
            path=Path("/root/samplesheet.csv"),
        ),
        "outdir": NextflowFileParameter(
            type=LatchOutputDir,
            default=LatchDir("latch://1721.account/funcscan-outputs"),
            path=Path("/root/funcscan-outputs"),
        ),
        # "amp_skip_amplify": NextflowParameter(type=bool, default=False),
        # "amp_skip_macrel": NextflowParameter(type=bool, default=False),
        # "amp_skip_ampir": NextflowParameter(type=bool, default=False),
        # "amp_ampir_model": NextflowParameter(type=str, default=""),
        # "amp_ampir_minlength": NextflowParameter(type=int, default=0),
        # "amp_skip_hmmsearch": NextflowParameter(type=bool, default=False),
        # "amp_hmmsearch_models": NextflowFileParameter(type=Optional[LatchDir], default=None),
        # "amp_hmmsearch_savealignments": NextflowParameter(type=bool, default=False),
        # "amp_hmmsearch_savetargets": NextflowParameter(type=bool, default=False),
        # "amp_hmmsearch_savedomains": NextflowParameter(type=bool, default=False),
        # "amp_ampcombi_db": NextflowFileParameter(type=Optional[LatchDir], default=None),
    },
)
