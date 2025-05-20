function struct_array = set_field_TD(struct_array, field_name, td_result)
names_split = strsplit(field_name,'.');
names_split = [names_split, 'LowerBack', 'TD'];
struct_array = setfield(struct_array,names_split{1:end},td_result);
% Construct the full field name dynamically
%full_field_name = strjoin([names_split, 'SU', 'LowerBack', 'TD'], '.');
%struct_array = setfield(struct_array, full_field_name, td_result);


