function [] = compareStructs(struct_A, struct_B, name_A, name_B)

arguments
    struct_A (1,1) struct
    struct_B (1,1) struct
    name_A (1,1) string = "struct_A"
    name_B (1,1) string = "struct_B"
end

fields_A = fieldnames(struct_A);
fields_B = fieldnames(struct_B);

fprintf("\n")

if ~isempty(setdiff(fields_A, fields_B)) || ~isempty(setdiff(fields_B, fields_A))
    fprintf("Additinal fields of %s are:\n", name_A)
    disp(setdiff(fields_A, fields_B))
    fprintf("Additinal fields of %s are:\n", name_B)
    disp(setdiff(fields_B, fields_A))
    fprintf("\n")
end

same_fields = fields_A(ismember(fields_A, fields_B));

for i = 1:length(same_fields)
    field = same_fields{i};
    if ~isequal(struct_A.(field), struct_B.(field))
        fprintf("Field %s is different.\n", field)
        fprintf("%s has value: %s\n", name_A, string(struct_A.(field)))
        fprintf("%s has value: %s\n", name_B, string(struct_B.(field)))
        fprintf("\n")
    end
end
