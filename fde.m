clear all
warning off
%%% QPSK over AWGN
%%%%% SIGNAL CONSTELLATION %%%%%%%%%
symbolBook=[1 -1];
bitBook=[0; 1];
n.BitPerSym=size(bitBook,2);
M=length(symbolBook);
n.SymPerFrame=1000;
h_l = [0.74 -0.514 0.37 0.216 0.062];
H_L = fft(h_l,n.SymPerFrame);
%%%%%%%%%%%%%% MONTE CARLO PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
n.BitsPerFrame=n.SymPerFrame*n.BitPerSym;
max_nFrame=2000;
fErrLim=100;
snr_db=0:20;
SYMBOLBOOK=repmat(transpose(symbolBook),1,n.SymPerFrame);
inputMatrix = load("fde_inputMatrix.mat");
inputMatrix = inputMatrix.inputMatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:height(inputMatrix)
    iv.cp = cell2mat(inputMatrix.cp(j));
    
    n.BitErrors=zeros(length(snr_db), 1);
    n.TransmittedFrames=zeros(length(snr_db), 1);
    n.ErroneusFrames=zeros(length(snr_db), 1);
    
    iv.chan_coefs = inputMatrix.f_l(j);
    iv.chan_coefs = cell2mat(iv.chan_coefs);
    iv.n1_plus_n2 = length(iv.chan_coefs)-1;
    
    
    for nEN = 1:length(snr_db) % SNR POINTS
        this_snr=snr_db(nEN);
        sigma_noise = 1/sqrt(10^(this_snr/10));
        if string(inputMatrix.equalizer(j)) == "mmse"
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
            G = G_transpose;
            w_head = inv(G*G'+(1/(10^(this_snr/10)))*I)*G*e;
        end
        while (n.TransmittedFrames(nEN)<max_nFrame) && (n.ErroneusFrames(nEN)<fErrLim)
            n.TransmittedFrames(nEN) = n.TransmittedFrames(nEN) + 1;
            %%%%%%%%%% INFORMATION GENERATION %%%%%%%%%%
            tr.SymIndices=randi(M,[1,n.SymPerFrame]);
            tr.SymVec=symbolBook(tr.SymIndices);
            tr.BitsMat=bitBook(tr.SymIndices,:)';
            
            if string(inputMatrix.equalizer(j)) == "mmse"
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
            else
                %%%%%%%%%%%%%ISI CHANNEL %%%%%%%%%%%%%%%%%%%%%
                tr.SymVec = [tr.SymVec(end-iv.cp+1:end) tr.SymVec];
                tr.SymVec_isi = conv(iv.chan_coefs,tr.SymVec); % filtering
                tr.SymVec_isi = tr.SymVec_isi(1:end-length(h_l)+1);
                %%%%%%%%%%%%%CHANNEL %%%%%%%%%%%%%%%%%%%%%
                noise=1/sqrt(2)*[randn(1, n.SymPerFrame+iv.cp) + 1j*randn(1,n.SymPerFrame+iv.cp)];
                recSigVec=tr.SymVec_isi+sigma_noise*noise;
                %%%%%%%%%%%%%EQUALIZATION %%%%%%%%%%%%%%%%%%%%%
                RECSIGVEC = fft(recSigVec(1+iv.cp:end));
                X_HEAD = RECSIGVEC./H_L;
                x_head = ifft(X_HEAD);
                recSigVec_equalized = x_head;
            end
            
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
        "displayname",strcat(...
        string(inputMatrix.cp(j)),"-cp",...
        " ",upper(string(inputMatrix.equalizer(j)))));
    %         ", ",string(iv.tap_no),"-taps",...
    hold on
end
legend("location","northeast");
axis square, grid on, grid minor;
set(gca,"FontSize",14);
ylabel("BER");
xlabel("SNR (dB)");