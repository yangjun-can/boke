%function [FIBFsLowToHigh FIBFsHighToLow EnergLekage]=FMD_Low2High_High2LowSacnning(xt,Fs,t)
%
% Calling sequence-
% [FIBFsLowToHigh FIBFsHighToLow EnergLekage]=FMD_Low2High_High2LowSacnning(xt,Fs,t)%
% All arguments are required
%
% Input-
%	  xt	      - 1-D time series
%	  Fs	      - Sampling Frequency
%	  t           - true time
%        
% Output-
%	  FIBFsLowToHigh      - FIBFs from Low to high frequency scan
%	  FIBFsHighToLow      - FIBFs from High to low frequency scan
%	  EnergLekage         - energy leakage
%     Two, Time-frequency-Energy plot
%
% References:   
%  Singh P, Joshi S.D., Patney R.K., Saha K., 
%  The Fourier Decomposition Method for nonlinear and nonstationary time
%  series analysis, arXiv:1503.06675 [stat.ME].


%
% code writer:  Pushpendra Singh, June, 2014.
% code writer:
%
% This is highly non optimzed initial vesrion of FDM, Matlab code. It can be optimzed.

function [xt_recov_IMFsLowToHigh, xt_recov_IMFsHighToLow, EnergyLekage] = FMD_Low2High_High2LowSacnning(xt,Fs,t)
    global PS_PhaseUnwrap;
    global CentralDiff;
    CentralDiff=1; % central finite difference, e.g. for delta function  中心有限差分，例如，对于delta函数  
    PS_PhaseUnwrap=0;
    threshold=-0*(10)^(-1); % idealy should be zero.
    L=length(xt);
    if 1 % adding zero at last, can create discontinuity in signal. e.g. org emd example 
        % 将信号长度固定为偶数 
        if rem(L,2) == 1 % odd
                L=L+1; % make it even for faster FFT
                %xt=[xt 0]; %  can create discontinuity
                xt=[xt xt(end)]; % repeat last data
        end
    end
    NFFT=L;
    N=(NFFT);
    Xk=fft(xt,NFFT)/L;     
    % k=1, Xk, is real, (2 to N/2) are complex, (N/2+1) is real, (N/2+2 to N) are complex conjugate of (N/2 to 2)    
    nTotalHormonics=N/2;
    
    %if 1 % this is much better, seen by examples 1.
        %%
        xt_recov_IMFsLowToHigh=zeros(1,N)';
        IFOfSignalIF=zeros(1,N)'; 
        xt_AnalyticFIBF=zeros(1,N)';
        xt_recov_IMFsLowToHigh(:,1)=(ones(1,N))*Xk(1); % first IMF DC
        init=2;
        p=2;     
        while(init<=N/2)            
            [kk xt_recov_FIBF xt_recov_AnalyticFIBF IFOfSignal]=getIMFsScanAllLowToHigh(Xk,Fs,init,nTotalHormonics,threshold);
            %% 
            init=kk+1
            xt_recov_IMFsLowToHigh(:,p)=xt_recov_FIBF;
            IFOfSignalIF(:,p-1)=IFOfSignal;
            xt_AnalyticFIBF(:,p-1)=xt_recov_AnalyticFIBF;            
            p=p+1;  
        end
         
        %subplot(2,1,1)
        sp_PlotTF(xt_AnalyticFIBF,IFOfSignalIF,t,Fs,Fs/2);
        title('FDM: Time-Frequncy-Energy estimate of FIBFs (LTH-FS)','FontSize',16,'FontName','Times');
        nn=(1:N);
        tttmp=Xk((N/2)+1).*cos(pi*(nn-1)); % N/2+1 part of FFT, last part of DFT
        xt_recov_IMFsLowToHigh(:,end+1)=tttmp'; % N/2+1 part of FFT, last part of DFT
    %end   
    
    
    %if 1 % this is much better, seen by examples 1 and 3.
        xt_recov_IMFsHighToLow=zeros(1,N)'; 
        IFOfSignalIF=zeros(1,N)'; 
        xt_AnalyticFIBF=zeros(1,N)';
        mm=1;        
        %tttmp=Xk((N/2)+1).*cos(pi*(1:N-1));
        xt_recov_IMFsHighToLow(:,mm)=tttmp'; % N/2+1 part of FFT, last part of DFT        
        final=N/2;        
        while(final>=2)
            mm=mm+1; 
            
            [kk xt_recov_FIBF xt_recov_AnalyticFIBF IFOfSignal]=getIMFsScanAllHighToLow(Xk,Fs,final,threshold);
            %% 
            final=kk-1
            xt_recov_IMFsHighToLow(:,mm)=xt_recov_FIBF;
            xt_AnalyticFIBF(:,mm-1)=xt_recov_AnalyticFIBF;
            IFOfSignalIF(:,mm-1)=IFOfSignal;            
        end
        %subplot(2,1,2)
        %figure 
        sp_PlotTF(xt_AnalyticFIBF,IFOfSignalIF,t,Fs,Fs/2);   
        title('FDM: Time-Frequncy-Energy estimate of FIBFs (HTL-FS)','FontSize',16,'FontName','Times');        
        %% residue
        xt_recov_IMFsHighToLow(:,mm+1)=(ones(1,N))*Xk(1); % first IMF DC        
    %end     
    
    %Calculate Energy Leakge
    FMF_energy=0;
    for np=1:length(xt_recov_IMFsHighToLow(1,:))         
         FMF_energy=FMF_energy+sum(xt_recov_IMFsHighToLow(:,np).*xt_recov_IMFsHighToLow(:,np));         
    end   
    SignalEnergy=sum(xt.*xt); 
    EnergyLekage=SignalEnergy-FMF_energy;          
    tt=1; 
    %%%  
    
    
    function [kk xt_recov_FIBF xt_recov_AnalyticFIBF IFOfSignal]=getIMFsScanAllLowToHigh(Xk,Fs,init,final,threshold)
    global PS_PhaseUnwrap;
    global CentralDiff;
    xt_recov_FIBF=0; 
    xt_recov_AnalyticFIBF=0;    
    N=length(Xk);
    %EndPoints=fix(length(Xk)*1/100);
    EndPoints=3; %threshold=-(10)^(-3);
    NbOfNegIF=(length(Xk)*1/100); % 1% is ok    
    Analytic_zt=0;
    n=1:length(Xk);
    col=1;
    PositiveIMF=1;
    for kk=init:final
        
         Analytic_zt=Analytic_zt+(Xk(kk)).*exp((1i*2*pi*(kk-1)*(n-1)/N)); 
         if PS_PhaseUnwrap
            IMFs_phase=sp_unwrap(angle(Analytic_zt));
         else
             IMFs_phase=unwrap(angle(Analytic_zt));
         end
         tmp1=(IMFs_phase);
         if CentralDiff
            %if 1 % central diff is better
            tmp=((tmp1(3:end)-tmp1(1:end-2))/2)*(Fs/(2*pi)); 
            tmp=[tmp(1) tmp tmp(end)];
            %else
            % five point diff
            %tmp=((8*tmp1(4:end-1)-8*tmp1(2:end-3)+tmp1(1:end-4)-tmp1(5:end))/12)*(Fs/(2*pi));
            %tmp=[tmp(1) tmp(1) tmp tmp(end) tmp(end)];
            %three point
            %tmp=((-3*tmp1(1:end-2)+4*tmp1(2:end-1)-tmp1(3:end))/2)*(Fs/(2*pi));
            %tmp=[tmp tmp(end) tmp(end)];
            %end
         else
            tmp=diff(tmp1)*(Fs/(2*pi)); % better derivative, e.g. use a wtooth & square waves as example
            tmp=[tmp tmp(end)];
            %tmp=[tmp(2) tmp(2:end-1) tmp(end-1) tmp(end-1)]; % First & last value may be negative, so 2:end-1
         end
         
         %if((sum(tmp(EndPoints:end-EndPoints)<0)>NbOfNegIF))
         if((min(tmp(EndPoints:end-(EndPoints-1)))<threshold) && (max(xt_recov_FIBF)~=0) && (PositiveIMF==1) && (kk<=final)) % check -ve IF
             PositiveIMF=0;             
             %break;
             recordIMF(:,col)=xt_recov_FIBF;
             recordAnalyticFIBFs(:,col)=xt_recov_AnalyticFIBF;
             recordIMF_IF(:,col)=xt_IF; % IF
             recordIMFkk(col)=kk-1;
             col=col+1;             
         end
         
         if(min(tmp(EndPoints:end-(EndPoints-1)))>=threshold) % check +ve IF, First value may be negative so 2:end
             PositiveIMF=1;             
         end         
         xt_recov_FIBF=2*real(Analytic_zt); % *2 as taking only half of component 
         xt_recov_AnalyticFIBF=2*Analytic_zt;
         xt_IF=tmp; % instantaneous Freq
         
         % First Value
         if(kk==init)
             recordIMF(:,col)=xt_recov_FIBF;
             recordAnalyticFIBFs(:,col)=xt_recov_AnalyticFIBF;
             recordIMF_IF(:,col)=xt_IF; % IF
             recordIMFkk(col)=kk;
             col=col+1;      
         end
         
         
         %% for last value
         if((kk==final) && (PositiveIMF==1))             
             recordIMF(:,col)=xt_recov_FIBF;
             recordAnalyticFIBFs(:,col)=xt_recov_AnalyticFIBF;
             recordIMF_IF(:,col)=xt_IF; % IF
             recordIMFkk(col)=kk;             
         end         
    end 
    getHighestValue=length(recordIMFkk);
    xt_recov_FIBF=recordIMF(:,getHighestValue); % return last value IMF
    xt_recov_AnalyticFIBF=recordAnalyticFIBFs(:,getHighestValue);
    IFOfSignal=recordIMF_IF(:,getHighestValue);
    kk=recordIMFkk(end); % return last value
    tt=1;
    
    function [kk xt_recov_FIBF xt_recov_AnalyticFIBF IFOfSignal]=getIMFsScanAllHighToLow(Xk,Fs,final,threshold)
    global PS_PhaseUnwrap;
    global CentralDiff;
    xt_recov_FIBF=0; 
    xt_recov_AnalyticFIBF=0;    
    N=length(Xk);
    %EndPoints=fix(length(Xk)*1/100);
    EndPoints=3; %threshold=-(10)^(-3);
    NbOfNegIF=(length(Xk)*1/100); % 1% is ok    
    Analytic_zt=0;
    n=1:length(Xk);
    col=1;
    PositiveIMF=1;
    for kk=final:-1:2
        
        Analytic_zt=Analytic_zt+(Xk(kk)).*exp((1i*2*pi*(kk-1)*(n-1)/N));        
        if PS_PhaseUnwrap
            IMFs_phase=sp_unwrap(angle(Analytic_zt));
        else
            IMFs_phase=unwrap(angle(Analytic_zt));
        end
         tmp1=(IMFs_phase);
         if CentralDiff
            tmp=((tmp1(3:end)-tmp1(1:end-2))/2)*(Fs/(2*pi));
            tmp=[tmp(1) tmp tmp(end)];
         else
            tmp=diff(tmp1)*(Fs/(2*pi)); % better derivative, e.g. use a wtooth & square waves as example
            tmp=[tmp tmp(end)];
         end
         
         %if((sum(tmp(EndPoints:end-EndPoints)<0)>NbOfNegIF))
         if( (min(tmp(EndPoints:end-(EndPoints-1)))<threshold) && (max(xt_recov_FIBF)~=0) && (PositiveIMF==1) && (kk>=2) ) % check -ve IF
         %if((sum(tmp(5:end-5)<0)>0) && (max(xt_recov_FIBF)~=0) && (PositiveIMF==1) && (kk>=2)) % check -ve IF
             PositiveIMF=0;             
             %break;
             recordIMF(:,col)=xt_recov_FIBF;
             recordAnalyticFIBFs(:,col)=xt_recov_AnalyticFIBF;
             recordIMF_IF(:,col)=xt_IF; % IF
             recordIMFkk(col)=kk+1;
             col=col+1;             
         end
         
         if(min(tmp(EndPoints:end-(EndPoints-1)))>=threshold) % check +ve IF, First value & last value may be negative so 2:end
         %if(min(tmp(5:end-5))>0)
             PositiveIMF=1;             
         end         
         xt_recov_FIBF=2*real(Analytic_zt); % *2 as taking only half of component
         xt_recov_AnalyticFIBF=2*Analytic_zt;
         xt_IF=tmp; % IF
         % First Value
         if(kk==final)
             recordIMF(:,col)=xt_recov_FIBF;
             recordAnalyticFIBFs(:,col)=xt_recov_AnalyticFIBF;
             recordIMF_IF(:,col)=xt_IF; % IF
             recordIMFkk(col)=kk;
             col=col+1;      
         end
         
         %% for last value
         if((kk==2) && (PositiveIMF==1))             
             recordIMF(:,col)=xt_recov_FIBF;
             recordAnalyticFIBFs(:,col)=xt_recov_AnalyticFIBF;
             recordIMF_IF(:,col)=xt_IF; % IF
             recordIMFkk(col)=kk;             
         end         
    end 
    getHighestValue=length(recordIMFkk);
    xt_recov_FIBF=recordIMF(:,getHighestValue); % return last value IMF
    xt_recov_AnalyticFIBF=recordAnalyticFIBFs(:,getHighestValue);
    IFOfSignal=recordIMF_IF(:,getHighestValue);
    kk=recordIMFkk(end); % return last value
    tt=1;    
    
    function sp_PlotTF(xt_recov_FIBFs,xt_recov_IMFs_IF,t,Fs,fw1)        
    xt_recov_FIBFs_abs=abs(xt_recov_FIBFs);
    %% To plot Use code of RCADA    
    %[nt,tscale,fscale]=PlotTF_FFT(freq,amp,t0,t1,fres,tres,fw0,fw1,tw0,tw1);
    [nt,tscale,fscale]=PlotTF_FFT(xt_recov_IMFs_IF(:,1:end),xt_recov_FIBFs_abs(:,1:end),t,Fs,fw1); % magnitude value 
    %[nt,tscale,fscale]=PlotTF_FFT(xt_recov_IMFs_IF(:,1:end),xt_recov_IMFs_AnalyticAbs(:,1:end),t,Fs,fw1); % magnitude value
    if 0
    imagesc(tscale,fscale,nt);
    axis xy;
    surf(tscale,fscale,nt); % for 3d graph
    end    
    q=fspecial('gaussian',7,0.6);  
    nsu=filter2(q,nt);
    nsu=filter2(q,nsu);
    figure;
    %subplot(5,1,1);
    imagesc(tscale,fscale,nsu.^.5); 
    axis xy;
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    %title('FMD based Time Frequncy Plot');
    title('FDM: Time-Frequency plot of FIBFs','FontSize',16,'FontName','Times'); 
    %% Use end code of RCADA
    
    FntSize=16;
    set(gca,'FontSize',FntSize,'FontName','Times')
    h1 = get(gca, 'xlabel');
    set(h1,'FontSize',FntSize,'FontName','Times')
    h1 = get(gca, 'ylabel');
    set(h1,'FontSize',FntSize,'FontName','Times') 
    
    
    
    function [nt,tscale,fscale]=PlotTF_FFT(freq,amp,t,Fs,fw1)
    
    t0=0;t1=t(end);    
    multfactor=4;
    if(length(freq(:,1))>=100*multfactor) 
        fres=100*multfactor; tres=100*multfactor;
    else
        fres=length(freq(:,1));    tres=fres;
    end
    
    %fres=268;    tres=fres; % for earthquake data
    
    %fw0=0; fw1=max(max(freq));
    fw0=min(min(freq)); fw1=max(max(freq));
    if(fw0<0)          fw0=0;    end
    %fw1=Fs/2; % max frequency in plot
    tw0=t0;     tw1=t1;
    %  4.call nspplote.m to plot the figure of time-frequency spectrum 
    %----- Get the values to plot the spectrum
    lscale=0;
    [nt,tscale,fscale] = nspplote(freq,amp,t0,t1,fres,tres,fw0,fw1,tw0,tw1,lscale); 
    %for clear appearance,clear the low-freqency if you want 
    setclear=0;
    if setclear==1
        nt(1,:)=0;
    end
    
    if 1% marginal spectrum plot
        %marginal spectrum
        ms=sum(nt,2);
        hsp_fre1=fw1; %give frequency axis max value:
        %comparisoncoefficients
        hco=Fs/2/hsp_fre1;
        figure;
        plot(fscale(1:end),hco*ms(1:end));
        %loglog(fscale(1:end),hco*ms(1:end));
        xlabel('Frequency (Hz)')
        ylabel('Spectral Density')
        title('Marginal spectrum by FDM','FontSize',16,'FontName','Times');    
        FntSize=16;
        set(gca,'FontSize',FntSize,'FontName','Times')
        h1 = get(gca, 'xlabel');
        set(h1,'FontSize',FntSize,'FontName','Times')
        h1 = get(gca, 'ylabel');
        set(h1,'FontSize',FntSize,'FontName','Times')
    
    end
    if 1% E(t) plot, Energy fluctuation with time
        %marginal spectrum
        ms=sum(nt,1);
        hsp_fre1=fw1; %give frequency axis max value:
        %comparisoncoefficients
        hco=Fs/2/hsp_fre1;
        figure;
        plot(tscale,hco*ms(1:end));
        %loglog(fscale(1:end),hco*ms(1:end));
        xlabel('Time (s)')
        ylabel('Energy Density')
        title('E(t) by FDM','FontSize',16,'FontName','Times');    
        FntSize=16;
        set(gca,'FontSize',FntSize,'FontName','Times')
        h1 = get(gca, 'xlabel');
        set(h1,'FontSize',FntSize,'FontName','Times')
        h1 = get(gca, 'ylabel');
        set(h1,'FontSize',FntSize,'FontName','Times')
    
    end
    
    