function [recon, magphase] = SpectrumInversion (magspect, magphase, fftSize_smp, grainSize_smp, grainStep_smp, phaseRateScale)
            
temp1 = phaseRateScale/grainSize_smp;

        %   check through the magnitude spectrum to find peak
            for j=2:1:fftSize_smp/2-1
                if(magspect(j)>magspect(j-1) && magspect(j)>magspect(j+1))

%  Phase estimate using quadratic interpolation as per Julios O. Smith
%  http://www.dsprelated.com/dspbooks/sasp/Quadratic_Interpolation_Spectral_Peaks.html
%  Use quadratic interpolation to estimate the real peak position

            alpha=magspect(j-1);
            beta=magspect(j);
            gamma=magspect(j+1);
            denom=alpha-2*beta+gamma;


            % initialize p
            
            if(denom~=0)
            p=0.5*(alpha-gamma)/denom;
            else
            p=0;
            end

        % checking for peaks

        phaseRate=2*pi*(j-1+p)*temp1;    %adjusted phase rate
        magphase(j)= magphase(j) + grainStep_smp*phaseRate; %phase accumulator for this peak bin
        peakPhase=magphase(j);
    

% Apply_input simple phase locking around the peaks.
% The phase relationships of the bins around the peak were determined by
% some simple experiments inspired by: 
% Laroche/Dolson "About This Phasiness Business" (1997), which mentions a paper
% M.S. Puckette "Phase-locked vocoder" (1995).  
% http://msp.ucsd.edu/Publications/mohonk95.pdf

% According to Laroche/Dolson:
% "Puckette in [5] recognized that for a constant-frequency_input constant-amplitude sinusoid 
% the sy_inputnthesis phases around the maximum of the Fourier transform should exhibit +/- pi 
% alternations and proposed a very_input simple way_input to constrain them to do so".


% Do 0/pi phase shift thing around the peak

% If actual peak is to the right
        if (p>0)

        % First bin to right has pi shift
            bin=j+1;
            magphase(bin)=peakPhase+pi;
            bin=j-1;
      
            while(bin>1 && magspect(bin)<magspect(bin+1))

    % Bins to left have shift of pi
                
                magphase(bin)=peakPhase+pi;
                bin=bin-1;
            end

            bin=j+2;
             while(bin<(fftSize_smp/2) && magspect(bin)<magspect(bin-1))

    % Other bins (apart from the immediate bin to the right of the peak)
    % have zero shift
    
                magphase(bin)=peakPhase+0;
                bin=bin+1;

             end        
        end
     
     
% Peak is to the left
      
     if(p<0)

    % First bin to left has pi shift
            bin=j-1;
            magphase(bin)=peakPhase+pi;

            bin=j+1;

 % similarly, like before, instead of using bin< fftSize_smp/2-1 
 % I use bin<fftSize_smp/2
        while(bin<(fftSize_smp/2) && magspect(bin)<magspect(bin-1))
            
% Bins to right have shift of pi
                magphase(bin)=peakPhase+pi;
                bin=bin+1;
            
        end
         bin=j-2;
         while(bin>1 && magspect(bin)<magspect(bin+1))
             
% Other bins to left have zero shift
            magphase(bin)=peakPhase+0;
            bin=bin-1;
        end
     end
                end
            end  

       recon=magspect.*exp(1i*magphase);       %reconstruct with new phase
       recon=[recon(1:fftSize_smp/2); zeros(fftSize_smp/2,1)];

   % Using conj for conjugate symmetry- reconstruct other half of FFT
   
    for k=2:1:fftSize_smp/2
    recon(fftSize_smp+2-k)=conj(recon(k));
    recon(fftSize_smp+2-k)=conj(recon(k));
    end

    % IFFT of the result and subsequent windowing
    
        recon(1)=0;
        recon(fftSize_smp/2+1)=0;
    
end
            
