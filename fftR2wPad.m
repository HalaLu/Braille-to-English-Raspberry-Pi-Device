function out = fftR2wPad (signal, size)

    if numel(signal) > size
        error('input signal larger than FFT size')
    elseif numel(signal) == size
        out = fft(signal,size);
    else
        temp = [signal ; zeros(size-numel(signal),1)];
        out = fft(temp,size);
    end  

end