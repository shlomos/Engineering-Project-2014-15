function [  ] = compute_norm_on_data( filename, GT_filename )
% This function estimates the given CT and GT and normal dist. and computes
% its expectation and std

global info;

CT = load_untouch_nii(filename);
CT_img = double(CT.img); %rows X cols X slices
GT = load_untouch_nii(GT_filename);
GT_img = double(GT.img); %rows X cols X slices

CT_img_row = CT_img(:);
GT_img_row = GT_img(:);

if (length(CT_img_row) ~= length(GT_img_row))
    fprintf('ERROR! the lengths of the CT and GT are not the same in case: %s\n', filename);
    size(CT_img)
    size(GT_img)
    return;
end

only_GT_on_CT_row = CT_img_row(logical(GT_img_row));

[mu,sigma,muci,sigmaci] = normfit(only_GT_on_CT_row);

fprintf('file: %s.\n', filename); 
fprintf('mu is: %f.\n', mu);
fprintf('sigms is: %f.\n\n', sigma);


info = [info; {filename, mu, sigma}];

end

