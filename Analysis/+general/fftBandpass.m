% function XF = fftBandpass(X, FS, FS1, FP1, FP2, FS2)
%
% Bandpass filter for the signal X (time x trials). An acuasal fft
% algorithm is applied (i.e. no phase shift). The filter functions is
% constructed from a Hamming window.
%
% Fs : sampling frequency
%
% The passbands (Fp1 Fp2) and stop bands (Fs1 Fs2) are defined as
%                 -----------
%                /           \
%               /             \
%              /               \
%             /                 \
%   ----------                   -----------------
%           Fs1  Fp1       Fp2  Fs2
%
% If no output arguments are assigned the filter function H(f) and
% impulse response are plotted.
%
% NOTE: for long data traces the filter is very slow.
%
%------------------------------------------------------------------------
% Ole Jensen, Brain Resarch Unit, Low Temperature Laboratory,
% Helsinki University of Technology, 02015 HUT, Finland,
% Report bugs to ojensen@neuro.hut.fi
%------------------------------------------------------------------------

%    Copyright (C) 2000 by Ole Jensen
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You can find a copy of the GNU General Public License
%    along with this package (4DToolbox); if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
function xf = fftBandpass(x,Fs,Fs1,Fp1,Fp2,Fs2)
    if size(x,1) == 1
        x = x';
    end
    % Make x even
    Norig = size(x,1);
    if rem(Norig,2)
        x = [x' zeros(size(x,2),1)]';
    end

    % Normalize frequencies
    Ns1 = Fs1/(Fs/2);
    Ns2 = Fs2/(Fs/2);
    Np1 = Fp1/(Fs/2);
    Np2 = Fp2/(Fs/2);

    % Construct the filter function H(f)
    N = size(x,1);
    Nh = N/2;

    B = fir2(N-1,[0 Ns1 Np1 Np2 Ns2 1],[0 0 1 1 0 0]);
    H = abs(fft(B));  % Make zero-phase filter function
    IPR = real(ifft(H));
    if nargout == 0
        subplot(2,1,1)
        f = Fs*(0:Nh-1)/(N);
        plot(f,H(1:Nh));
        xlim([0 2*Fs2])
        ylim([0 1]);
        title('Filter function H(f)')
        xlabel('Frequency (Hz)')
        subplot(2,1,2)
        plot((1:Nh)/Fs,IPR(1:Nh))
        xlim([0 2/Fp1])
        xlabel('Time (sec)')
        ylim([min(IPR) max(IPR)])
        title('Impulse response')
    end


    if size(x,2) > 1
        for k=1:size(x,2)
            xf(:,k) = real(ifft(fft(x(:,k)) .* H'));
        end
        xf = xf(1:Norig,:);
    else
        xf = real(ifft(fft(x') .* H));
        xf = xf(1:Norig);
    end
end