#!/usr/bin/env sh

# jmodeltest
java -jar ~/Software/jmodeltest2/dist/jModelTest.jar -d pacus.haps.trimmed.nex -s 5 -S NNI -g 4 -tr 8 -f -AICc -BIC -w -o pacus.haps.trimmed.nex.out

# beast run1 ~1.42 minutes/million states
beast -beagle_auto -overwrite -seed 1408191 pacus.haps.trimmed.run1.xml

# beast run2
beast -beagle_auto -overwrite -seed 1408192 pacus.haps.trimmed.run2.xml

# beast run3
beast -beagle_auto -overwrite -seed 1408193 pacus.haps.trimmed.run3.xml
