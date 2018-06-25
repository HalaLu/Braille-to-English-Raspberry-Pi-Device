function [grain] = getgrain_seq(snd, grainSize_smp, filePos_smp, spray_smp, sprayOffset_smp, sequenceRate, sndLength_smp, idx)
    
    pos = filePos_smp + round(idx * sequenceRate * grainSize_smp) + round(rand * spray_smp) - sprayOffset_smp + 1;
    pos = rem(pos,(sndLength_smp - grainSize_smp));
    if pos < 0
        pos = sndLength_smp - grainSize_smp + pos;
    end
    grain = snd(pos:(pos + grainSize_smp - 1),:);

end