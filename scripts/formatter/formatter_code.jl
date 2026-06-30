using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
Pkg.update()

using JuliaFormatter

targets = isempty(ARGS) ? ["./src"] : ARGS

fmt_opts = (;
    whitespace_ops_in_indices = true,
    remove_extra_newlines = true,
    verbose = true,
    always_for_in = true,
    whitespace_typedefs = true,
    conditional_to_if = true,
    join_lines_based_on_source = true,
    separate_kwargs_with_semicolon = true,
    format_markdown = true,
    ignore = ["*LICENSE.md", "how_to/install.md"], # install has complicated formatting
    # always_use_return = true, # Disabled since it throws a lot of false positives
)

function format_target(target)
    if isfile(target)
        @show file_path = abspath(target)
        format(file_path; fmt_opts...)
    elseif isdir(target)
        for (root, dir, files) in walkdir(target)
            for f in files
                !((occursin(".jl", f) || occursin(".md", f))) && continue
                @show file_path = abspath(root, f)
                format(file_path; fmt_opts...)
            end
        end
    else
        @warn "Path not found: $target"
    end
end

for target in targets
    format_target(target)
end