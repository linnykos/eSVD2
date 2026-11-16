#!/bin/zsh
# Green split-half study: steps 3 and 4, started once the four fits of
# step 2 have ended. Drafted by Claude for Kevin Z. Lin, 2026-11-16.
# Run from the repository root; logs go to output/ beside this script.
O="$HOME/Library/CloudStorage/Dropbox/Collaboration-and-People/Kathryn Roeder - private/eSVD2/out/2026-11-16_green"
D=additional_context/2026-11-16_green
while pgrep -f "02_fit-esvd_claude.R" > /dev/null; do /bin/sleep 60; done
for l in split1_half1 split1_half2 split2_half1 split2_half2; do
  if [ -f "$O/green_${l}_esvd.RData" ]; then
    Rscript $D/03_apply-rules_claude.R $l > $D/output/rules_$l.log 2>&1
    echo "$(date) $l rules exit $?"
  else
    echo "$(date) $l: NO FIT"
  fi
done
Rscript $D/04_evaluate_claude.R > $D/output/evaluate.log 2>&1
echo "$(date) evaluate exit $?"
