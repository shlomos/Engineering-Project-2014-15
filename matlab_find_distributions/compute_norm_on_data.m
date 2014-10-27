function [  ] = compute_norm_on_data( filename )
% This function estimates the given CT and GT and normal dist. and computes
% its expectation and std

CT = load_untouch_nii('C:\Users\moshesamson\Desktop\cropped_liver\case9\case9-15-07-2007.nii');
CT_img = double(CT.img); %rows X cols X slices
GT = load_untouch_nii('C:\Users\moshesamson\Desktop\cropped_liver\case9\case9-15-07-2007_GT.nii');
GT_img = double(GT.img); %rows X cols X slices

CT_img_row = CT_img(:);
GT_img_row = GT_img(:);

only_GT_on_CT_row = CT_img_row(logical(GT_img_row));

[mu,sigma,muci,sigmaci] = normfit(only_GT_on_CT_row);

disp('file: ') , filename
disp('mu is: '), mu
disp('sigms is: '), sigma


end

