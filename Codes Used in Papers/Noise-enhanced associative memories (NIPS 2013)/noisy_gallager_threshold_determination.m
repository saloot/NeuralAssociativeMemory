n = 400;
m = 200;
d_v = 4;
d_c = 8;
upsilon = 0.0;
nu = 0.0;
varphi = 0.75;
itr_max = 100;
no_tested_instances = 100;
successful_itr = [];
BER_total = [];
decoding_thr = 0;
for p_e = .05:.1:1
    [BER,itr] = noisy_gallager(n,m,d_v,d_c,upsilon,nu,varphi,itr_max,p_e,no_tested_instances);

    BER_total = [BER_total,BER];
    successful_itr = [successful_itr,itr];
    if (BER/p_e > .8)
        decoding_thr = p_e;
        break;
    end
end

% fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Noisy_Gallagher/Decoding_thr_n_',num2str(n),'_m_',num2str(m),'_dv_',num2str(d_v),'_dc_',num2str(d_c),...
%     '_upsilon_',num2str(upsilon),'_nu_',num2str(nu),'_varphi_',num2str(varphi),'_itr_',num2str(itr_max),'.txt'],'w');
%     fprintf(fid, '%d %d', [2 job_id]);    
% fclose(fid);

save(['/scratch/amir/Fault_Tolerant_Decoding/Noisy_Gallagher/Decoding_thr_n_',num2str(n),'_m_',num2str(m),'_dv_',num2str(d_v),'_dc_',num2str(d_c),...
     '_upsilon_',num2str(upsilon),'_nu_',num2str(nu),'_varphi_',num2str(varphi),'_itr_',num2str(itr_max),'.mat'],'decoding_thr','successful_itr');