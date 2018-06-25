function [grain] = getgrain_rnd(snd, grainSize_smp)
    
    pos = round(rand * (size(snd,1) - grainSize_smp - 1)) + 1;
    grain = snd(pos:(pos + grainSize_smp - 1),:);

end



