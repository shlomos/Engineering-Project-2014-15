function [ ] = function1( )
% This function receives as input from the user the location of the CT and
% GT files, and compute the normal distribution values for each example.
global info;
info = {'filename', 'mu', 'sigma'};

main_dir = input('Please insert the path of the main directory: \n', 's');
list_of_examples = dir(main_dir);

 for i = 3:size(list_of_examples,1)
     
     % there are no GT for cases 11-17
     if ((strcmp(list_of_examples(i).name, 'case11'))|| ...
         (strcmp(list_of_examples(i).name, 'case12'))|| ...
         (strcmp(list_of_examples(i).name, 'case13'))|| ...
         (strcmp(list_of_examples(i).name, 'case14'))|| ...
         (strcmp(list_of_examples(i).name, 'case15'))|| ...
         (strcmp(list_of_examples(i).name, 'case16'))|| ...
         (strcmp(list_of_examples(i).name, 'case17')))
         continue
     end
     
     list_files_of_case = dir(strcat(main_dir , '\', list_of_examples(i).name));
     for j = 3:2:size(list_files_of_case,1)
         filename = strcat(main_dir, '\', list_of_examples(i).name, '\', list_files_of_case(j).name);
         GT_filename = strcat(main_dir, '\', list_of_examples(i).name, '\', list_files_of_case(j+1).name);
         compute_norm_on_data(filename, GT_filename, 1);
     end
 end
 
 fprintf('The mean of mu is   : %f\n',mean(cell2mat(info(2:end, 2))))
 fprintf('The mean of sigma is: %f\n',mean(cell2mat(info(2:end, 3))))
 xlswrite('output.xls', info);

end

