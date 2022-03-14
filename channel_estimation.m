clc, clear;
warning off
girdi = 0;
girmedi = 0;
%% SIGNAL CONSTELLATION %%%%%%%%%
symbolBook=[1 -1];
bitBook=[0; 1];
n.BitPerSym=size(bitBook,2);
M=length(symbolBook);

h = [0.74 -0.514 0.37 0.216 0.062];
iv.chan_coefs = h;
L = length(h);
n1 = 0;
h_head_lms = [Inf Inf Inf];

%%%%%%%%%%%%%% MONTE CARLO PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%

max_nFrame=2000;
fErrLim=100;
snr_db=0:2:30;

inputMatrix = load("fde_inputMatrix.mat");
inputMatrix = inputMatrix.inputMatrix;

for j = 1:height(inputMatrix)-1
    no_of_pilots = inputMatrix.pilots(j);
    n2 = no_of_pilots;
    n.SymPerFrame = 1000 + no_of_pilots;
    n.BitsPerFrame=n.SymPerFrame*n.BitPerSym;
    SYMBOLBOOK=repmat(transpose(symbolBook),1,n.SymPerFrame-no_of_pilots);
    %%
    n.BitErrors=zeros(length(snr_db), 1);
    n.TransmittedFrames=zeros(length(snr_db), 1);
    n.ErroneusFrames=zeros(length(snr_db), 1);
    
    for nEN = 1:length(snr_db) % SNR POINTS
        this_snr=snr_db(nEN);
        sigma_noise = 1/sqrt(10^(this_snr/10));
        while (n.TransmittedFrames(nEN)<max_nFrame) && (n.ErroneusFrames(nEN)<fErrLim)
            n.TransmittedFrames(nEN) = n.TransmittedFrames(nEN) + 1;
            %% INFORMATION GENERATION %%%%%%%%%%
            tr.SymIndices=randi(M,[1,n.SymPerFrame]);
            tr.SymVec=symbolBook(tr.SymIndices);
            tr.BitsMat=bitBook(tr.SymIndices,:)';
            
            %% ISI CHANNEL %%%%%%%%%%%%%%%%%%%%%
            tr.SymVec_isi = conv(h,tr.SymVec); % filtering
            tr.SymVec_isi = tr.SymVec_isi(1:length(tr.SymVec_isi)-L+1);
            %% CHANNEL %%%%%%%%%%%%%%%%%%%%%
            noise=1/sqrt(2)*[randn(1, n.SymPerFrame) + 1j*randn(1,n.SymPerFrame)];
            recSigVec=tr.SymVec+sigma_noise*noise;
            %%
            R_z = zeros(L,L);
            s = tr.SymVec;
            for jj = 1:L
                for k = 1:L
                    for l = 5:n2
                        R_z(jj,k) = R_z(jj,k) + s(l-k+1)*conj(s(l-jj+1));
                    end
                end
            end
            
            r_sz = zeros(L,1);
            z = recSigVec;
            for jj = 1:L
                for k = n1+L:n2
                    r_sz(jj) = r_sz(jj) + z(k).*conj(s(k-jj+1));
                end
            end
            h_head = inv(R_z)*real(r_sz);
            h_head = h_head';
            %%
            if (cond(R_z) < 10^-4) || (cond(R_z) > 10^4)
                girdi = girdi + 1;
                n.TransmittedFrames(nEN) = n.TransmittedFrames(nEN) - 1;
                if no_of_pilots == 5
                    no_of_pilots = no_of_pilots*2;
                    n2 = no_of_pilots;
                    n.SymPerFrame = 2000 + no_of_pilots;
                    n.BitsPerFrame=n.SymPerFrame*n.BitPerSym;
                    SYMBOLBOOK=repmat(transpose(symbolBook),1,n.SymPerFrame-no_of_pilots);
                end
                continue;
            end
            girmedi = girmedi + 1;
            h_head_lms(j) = min([h_head_lms(j) sqrt(mean((h_head - h).^2))]);
            %%
            iv.n1_plus_n2 = length(h_head)-1;
            iv.tap_no = inputMatrix.tap_no(j);
            iv.m1_plus_m2 = iv.tap_no - 1;
            e = zeros(iv.n1_plus_n2 + iv.m1_plus_m2 + 1, 1);
            e(1) = 1;
            G = zeros(iv.n1_plus_n2 + iv.m1_plus_m2 + 1, iv.m1_plus_m2 + 1);
            for i=1:iv.m1_plus_m2+1
                G(:,i) = [ zeros(i-1,1);h_head';zeros(iv.m1_plus_m2-i+1,1)];
            end
            G_transpose = G';
            I = eye(iv.m1_plus_m2 + 1);
            switch string(inputMatrix.equalizer(j))
                case "mmse"
                    G = G_transpose;
                    w_head = inv(G*G'+(1/(10^(this_snr/10)))*I)*G*e;
            end
            %%
            recSigVec = recSigVec(no_of_pilots+1:end);
            tr.BitsMat = tr.BitsMat((no_of_pilots+1):end);
            %% EQUALIZATION %%%%%%%%%%%%%%%%%%%%%
            recSigVec_equalized = conv(w_head,recSigVec);
            recSigVec_equalized = recSigVec_equalized(1:length(recSigVec_equalized)...
                -(iv.tap_no-1));
            %% DETECTOR %%%%%%%%%%%%
            RECSIGVEC=repmat(recSigVec_equalized,length(symbolBook),1);
            distance_mat=abs(SYMBOLBOOK-RECSIGVEC);
            [~, det_sym_ind]=min(distance_mat,[],1);
            detected_bits=[bitBook(det_sym_ind, :)]';
            err = sum(sum(abs(tr.BitsMat-detected_bits)));
            n.BitErrors(nEN)=n.BitErrors(nEN)+err;
            if err~=0
                n.ErroneusFrames(nEN)=n.ErroneusFrames(nEN)+1;
            end
        end % End of while loop
        sim_res=[n.BitErrors n.TransmittedFrames]
    end %end for (SNR points)
    inputMatrix.sim_res(j) = struct("nBitErrors",n.BitErrors,"nTransmittedFrames",n.TransmittedFrames);
    semilogy(snr_db, n.BitErrors./n.TransmittedFrames/n.BitsPerFrame,...
        strcat(string(inputMatrix.color(j)),...
        string(inputMatrix.linestyle(j)),...
        string(inputMatrix.marker(j))),...
        "linewidth", 2,...
        "displayname",strcat(...
        string(inputMatrix.pilots(j)),"-pilot",...
        ", ",string(iv.tap_no),"-taps",...
        ", ",upper(string(inputMatrix.equalizer(j)))));
    hold on
end
%%
clear all
warning off
%%%%% SIGNAL CONSTELLATION %%%%%%%%%
symbolBook=[1 -1];
bitBook=[0; 1];
n.BitPerSym=size(bitBook,2);
M=length(symbolBook);

f_l_1 = [0.74 -0.514 0.37 0.216 0.062];
%%%%%%%%%%%%%% MONTE CARLO PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
n.SymPerFrame=3000;
n.BitsPerFrame=n.SymPerFrame*n.BitPerSym;
max_nFrame=2000;
fErrLim=100;
snr_db=0:20;
SYMBOLBOOK=repmat(transpose(symbolBook),1,n.SymPerFrame);
inputMatrix = load("inputMatrix.mat");
inputMatrix = inputMatrix.inputMatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 4
    
    n.BitErrors=zeros(length(snr_db), 1);
    n.TransmittedFrames=zeros(length(snr_db), 1);
    n.ErroneusFrames=zeros(length(snr_db), 1);
    
    iv.chan_coefs = inputMatrix.f_l(j);
    iv.chan_coefs = cell2mat(iv.chan_coefs);
    iv.n1_plus_n2 = length(iv.chan_coefs)-1;
    iv.tap_no = inputMatrix.tap_no(j);
    iv.m1_plus_m2 = iv.tap_no - 1;
    e = zeros(iv.n1_plus_n2 + iv.m1_plus_m2 + 1, 1);
    e(1) = 1;
    G = zeros(iv.n1_plus_n2 + iv.m1_plus_m2 + 1, iv.m1_plus_m2 + 1);
    for i=1:iv.m1_plus_m2+1
        G(:,i) = [ zeros(i-1,1);iv.chan_coefs';zeros(iv.m1_plus_m2-i+1,1)];
    end
    G_transpose = G';
    I = eye(iv.m1_plus_m2 + 1);
    for nEN = 1:length(snr_db) % SNR POINTS
        this_snr=snr_db(nEN);
        sigma_noise = 1/sqrt(10^(this_snr/10));
        switch string(inputMatrix.equalizer(j))
            case "zfe"
                w_head = inv(G'*G)*G'*e;
            case "mmse"
                G = G_transpose;
                w_head = inv(G*G'+(1/(10^(this_snr/10)))*I)*G*e;
        end
        
        while (n.TransmittedFrames(nEN)<max_nFrame) && (n.ErroneusFrames(nEN)<fErrLim)
            n.TransmittedFrames(nEN) = n.TransmittedFrames(nEN) + 1;
            %%%%%%%%%% INFORMATION GENERATION %%%%%%%%%%
            tr.SymIndices=randi(M,[1,n.SymPerFrame]);
            tr.SymVec=symbolBook(tr.SymIndices);
            tr.BitsMat=bitBook(tr.SymIndices,:)';
            
            %%%%%%%%%%%%%ISI CHANNEL %%%%%%%%%%%%%%%%%%%%%
            tr.SymVec_isi = conv(iv.chan_coefs,tr.SymVec); % filtering
            tr.SymVec_isi = tr.SymVec_isi(1:length(tr.SymVec_isi)...
                -(length(iv.chan_coefs)-1));
            %%%%%%%%%%%%%CHANNEL %%%%%%%%%%%%%%%%%%%%%
            noise=1/sqrt(2)*[randn(1, n.SymPerFrame) + 1j*randn(1,n.SymPerFrame)];
            recSigVec=tr.SymVec_isi+sigma_noise*noise;
            %%%%%%%%%%%%%EQUALIZATION %%%%%%%%%%%%%%%%%%%%%
            recSigVec_equalized = conv(w_head,recSigVec);
            recSigVec_equalized = recSigVec_equalized(1:length(recSigVec_equalized)...
                -(iv.tap_no-1));
            %%%% DETECTOR %%%%%%%%%%%%
            RECSIGVEC=repmat(recSigVec_equalized,length(symbolBook),1);
            distance_mat=abs(SYMBOLBOOK-RECSIGVEC);
            [~, det_sym_ind]=min(distance_mat,[],1);
            detected_bits=[bitBook(det_sym_ind, :)]';
            err = sum(sum(abs(tr.BitsMat-detected_bits)));
            n.BitErrors(nEN)=n.BitErrors(nEN)+err;
            if err~=0
                n.ErroneusFrames(nEN)=n.ErroneusFrames(nEN)+1;
            end
        end % End of while loop
        sim_res=[n.BitErrors n.TransmittedFrames]
    end %end for (SNR points)
    inputMatrix.sim_res(j) = struct("nBitErrors",n.BitErrors,"nTransmittedFrames",n.TransmittedFrames);
    semilogy(snr_db, n.BitErrors./n.TransmittedFrames/n.BitsPerFrame,...
        strcat(string(inputMatrix.color(j)),...
        string(inputMatrix.linestyle(j)),...
        string(inputMatrix.marker(j))),...
        "linewidth", 2,...
        "displayname","0-pilot, 10-taps, MMSE");
    hold on
end
legend("location","northeast");
axis square, grid on, grid minor,
set(gca,"FontSize",14);
ylabel("BER");
xlabel("SNR (dB)");