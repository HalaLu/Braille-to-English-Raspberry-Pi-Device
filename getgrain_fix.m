function [grain] = getgrain_fix(snd, grainSize_smp, filePos_smp, spray_smp, sprayOffset_smp)

    pos = filePos_smp + (round(rand * spray_smp)) + 1 - sprayOffset_smp;
    grain = snd(pos:(pos + grainSize_smp - 1),:);

end