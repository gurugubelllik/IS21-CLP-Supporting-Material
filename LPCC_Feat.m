function feat=LPCC_Feat(wav,fs,lporder,pre_emp_coeff)

[~,cep]=LPCC(wav,fs,lporder,pre_emp_coeff);

% del_1=deltas(cep',5);

% del_2=deltas(del_1,5);

y1=cep;% del_1' del_2'];

feat=cmvn(y1')';
end
