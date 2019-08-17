#!/usr/bin/env sh

# jmodeltest (run from temp-local-only)
java -jar ~/Software/jmodeltest2/dist/jModelTest.jar -d pacus.haps.trimmed.nex -s 5 -S NNI -g 4 -tr 8 -f -AICc -BIC -w -o pacus.haps.trimmed.nex.out

# beast run1 (run from temp-local-only)
beast -beagle_auto -overwrite -seed 1608191 pacus.haps.trimmed.run1.xml

# beast run2 (run from temp-local-only)
beast -beagle_auto -overwrite -seed 1608192 pacus.haps.trimmed.run2.xml

# beast run3 (run from temp-local-only)
beast -beagle_auto -overwrite -seed 1608193 pacus.haps.trimmed.run3.xml

# check with tracer (run from temp-local-only)
tracer pacus.haps.trimmed.run1.log pacus.haps.trimmed.run2.log pacus.haps.trimmed.run3.log

# combine the trees (run from scripts)
./tree-combiner.sh ../temp-local-only/pacus.haps.trimmed.run1.trees ../temp-local-only/pacus.haps.trimmed.run2.trees ../temp-local-only/pacus.haps.trimmed.run3.trees ../temp-local-only/pacus.haps.trimmed.combined.trees 336 ../temp-local-only/pacus.haps.trimmed.combined.tre

# check with figtree (run from scripts)
figtree ../temp-local-only/pacus.haps.trimmed.combined.tre

# copy to data dir
cp ../temp-local-only/pacus.haps.trimmed.combined.tre ../data/pacus.COI.tre
