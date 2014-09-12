 function [x,mc,mn,mx]=melbankm(p,n,fs,fl,fh,w)
% 0002 %MELBANKM determine matrix for a mel-spaced filterbank [X,MN,MX]=(P,N,FS,FL,FH,W)
% 0003 %
% 0004 % Inputs:    p   number of filters in filterbank or the filter spacing in k-mel [default 0.06]
% 0005 %        n   length of fft
% 0006 %        fs  sample rate in Hz
% 0007 %        fl  low end of the lowest filter as a fraction of fs (default = 0)
% 0008 %        fh  high end of highest filter as a fraction of fs (default = 0.5)
% 0009 %        w   any sensible combination of the following:
% 0010 %              't'  triangular shaped filters in mel domain (default)
% 0011 %              'n'  hanning shaped filters in mel domain
% 0012 %              'm'  hamming shaped filters in mel domain
% 0013 %
% 0014 %              'z'  highest and lowest filters taper down to zero (default)
% 0015 %              'y'  lowest filter remains at 1 down to 0 frequency and
% 0016 %               highest filter remains at 1 up to nyquist freqency
% 0017 %
% 0018 %               If 'ty' or 'ny' is specified, the total power in the fft is preserved.
% 0019 %
% 0020 % Outputs:    x     a sparse matrix containing the filterbank amplitudes
% 0021 %                  If the mn and mx outputs are given then size(x)=[p,mx-mn+1]
% 0022 %                 otherwise size(x)=[p,1+floor(n/2)]
% 0023 %                 Note that teh peak filter values equal 2 to account for the power
% 0024 %                 in the negative FFT frequencies.
% 0025 %           mc    the filterbank centre frequencies in mel
% 0026 %            mn    the lowest fft bin with a non-zero coefficient
% 0027 %            mx    the highest fft bin with a non-zero coefficient
% 0028 %                 Note: you must specify both or neither of mn and mx.
% 0029 %
% 0030 % Usage:    f=fft(s);            f=fft(s);
% 0031 %        x=melbankm(p,n,fs);        [x,mc,na,nb]=melbankm(p,n,fs);
% 0032 %        n2=1+floor(n/2);        z=log(x*(f(na:nb)).*conj(f(na:nb)));
% 0033 %        z=log(x*abs(f(1:n2)).^2);
% 0034 %        c=dct(z); c(1)=[];
% 0035 %
% 0036 % To plot filterbanks e.g.    n=256; fs=8000; plot((0:floor(n/2))*fs/n,melbankm(20,n,fs)')
% 0037 %
% 0038 % [1] S. S. Stevens, J. Volkman, and E. B. Newman. A scale for the measurement
% 0039 %     of the psychological magnitude of pitch. J. Acoust Soc Amer, 8: 185�19, 1937.
% 0040 % [2] S. Davis and P. Mermelstein. Comparison of parametric representations for
% 0041 %     monosyllabic word recognition in continuously spoken sentences.
% 0042 %     IEEE Trans Acoustics Speech and Signal Processing, 28 (4): 357�366, Aug. 1980.
% 0043 
% 0044 
% 0045 %      Copyright (C) Mike Brookes 1997
% 0046 %      Version: $Id: melbankm.m,v 1.7 2009/10/19 10:19:40 dmb Exp $
% 0047 %
% 0048 %   VOICEBOX is a MATLAB toolbox for speech processing.
% 0049 %   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
% 0050 %
% 0051 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0052 %   This program is free software; you can redistribute it and/or modify
% 0053 %   it under the terms of the GNU General Public License as published by
% 0054 %   the Free Software Foundation; either version 2 of the License, or
% 0055 %   (at your option) any later version.
% 0056 %
% 0057 %   This program is distributed in the hope that it will be useful,
% 0058 %   but WITHOUT ANY WARRANTY; without even the implied warranty of
% 0059 %   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% 0060 %   GNU General Public License for more details.
% 0061 %
% 0062 %   You can obtain a copy of the GNU General Public License from
% 0063 %   http://www.gnu.org/copyleft/gpl.html or by writing to
% 0064 %   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
% 0065 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 if nargin < 6
     w='tz'; % default options
     if nargin < 5
         fh=0.5; % max freq is the nyquist
         if nargin < 4
             fl=0; % min freq is DC
         end
     end
 end
 f0=700/fs;
 fn2=floor(n/2);     % bin index of Nyquist term
 if isempty(p)
     p=0.06;         % spacing = 0.06 kmel
 end
 if p<1
     p=round(log((f0+fh)/(f0+fl))/(p*log(17/7)))-1;
 end
 lr=log((f0+fh)/(f0+fl))/(p+1);
 % convert filter edges to fft bin numbers (0 = DC)
 bl=n*((f0+fl)*exp([0 1 p p+1]*lr)-f0);  % bins: [filter1-low filter1-mid filterp-mid filterp-high]
 b2=ceil(bl(2));
 b3=floor(bl(3));
 mc=(log(fl/f0+1)+(1:p)*lr)*1000/log(17/7);          % mel centre frequencies
 if any(w=='y')          % preserve power in FFT
     pf=log((f0+(b2:b3)/n)/(f0+fl))/lr;
     fp=floor(pf);
     r=[ones(1,b2) fp fp+1 p*ones(1,fn2-b3)];
     c=[1:b3+1 b2+1:fn2+1];
     v=2*[0.5 ones(1,b2-1) 1-pf+fp pf-fp ones(1,fn2-b3-1) 0.5];
     mn=1;
     mx=fn2+1;
 else
     b1=floor(bl(1))+1;            % lowest FFT bin required (0 = DC)
     b4=min(fn2,ceil(bl(4)))-1;    % highest FFT bin required (0 = DC)
     pf=log((f0+(b1:b4)/n)/(f0+fl))/lr;  % maps FFT bins to filter
     if pf(end)>p
         pf(end)=[];
         b4=b4-1;
     end
     fp=floor(pf);                  % FFT bin i contributes to filters fp(1+i-b1)+[0 1]
     pm=pf-fp;
     k2=b2-b1+1;
     k3=b3-b1+1;
     k4=b4-b1+1;
     r=[fp(k2:k4) 1+fp(1:k3)];
     c=[k2:k4 1:k3];
     v=2*[1-pm(k2:k4) pm(1:k3)];
     mn=b1+1;
     mx=b4+1;
 end
 if any(w=='n')
     v=1-cos(v*pi/2);      % convert triangles to Hanning
 elseif any(w=='m')
     v=1-0.92/1.08*cos(v*pi/2);  % convert triangles to Hamming
 end
 if nargout > 2
     x=sparse(r,c,v);
     if nargout == 3
         mc=mn;    % delete mc output for legacy code compatibility
         mn=mx;
     end
 else
     x=sparse(r,c+mn-1,v,p,1+fn2);
 end
 if ~nargout
     ng=3;
     if any(w=='n') || any(w=='m')
         ng=201;
     end
     fe=((f0+fl)*exp((0:p+1)*lr)-f0)'*fs;
     x=repmat(linspace(0,1,ng),p,1).*repmat(fe(3:end)-fe(1:end-2),1,ng)+repmat(fe(1:end-2),1,ng);
     v=2-abs(linspace(-2,2,ng));
     if any(w=='n')
         v=1-cos(v*pi/2);      % convert triangles to Hanning
    elseif any(w=='m')
         v=1-0.92/1.08*cos(v*pi/2);  % convert triangles to Hamming
    end
    v=repmat(v,p,1);
     if any(w=='y')
         v(1,1:floor(ng/2))=2;
         v(end,ceil(ng/2):end)=2;
     end
     plot(x'/1000,v','b');
     set(gca,'xlim',[fe(1) fe(end)]/1000);
     xlabel('Frequency (kHz)');
 end