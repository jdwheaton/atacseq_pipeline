#!/bin/bash

singularity exec -H $PWD -B /work/mc394/ docker://ewels/multiqc multiqc .