using Printf

function replace_template(template::String, target::String, replacements::Dict{String,String})
    """
    Create a target file by replacing keys with values from `replacements` in the template.

    Arguments:
        template::String: Path to the template file.
        target::String: Path to the output file.
        replacements::Dict{String,String}: Pairs of {old => new} substrings.
    """
    open(template, "r") do fin
        open(target, "w") do fout
            for line in eachline(fin)
                for (old, new) in replacements
                    line = replace(line, old => new)
                end
                println(fout, line)
            end
        end
    end
end

# Helper to convert list of strings to floats
function to_float(ls)
    return [parse(Float64, s) for s in ls]
end

# Simple argument parser using only ARGS
function parse_args()
    filename = "init.inp"
    scale = 1.0
    zscale_flag = false

    i = 1
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "-F" || arg == "--file"
            i += 1
            filename = ARGS[i]
        elseif arg == "-L" || arg == "--scale"
            i += 1
            scale = parse(Float64, ARGS[i])
        elseif arg == "--zscale"
            zscale_flag = true
        else
            println("Unknown argument: $arg")
        end
        i += 1
    end

    zscale = zscale_flag ? scale : 1.0

    return filename, scale, zscale
end

# Get arguments
filename, scale, zscale = parse_args()

println("Scaling file $filename with scale $scale and Z-scale $zscale")

replace_template("template.geo", "current.geo", Dict("NX"=>string(scale)))

# make the mesh file
run(`gmsh current.geo -order 1 -3 -o init.inp`)

# Open mesh file and read lines
lines = readlines(filename)
println(lines[1])

# Find where the element definition starts
elem_def = 0
for (k, line) in enumerate(lines)
    if occursin("*ELEMENT, type=C3D4,", line)
        println("Found solid element definition at line $k: $line")
        global elem_def = k
        break
    end
end

eof = length(lines)
println("Element definition starts at line $elem_def and ends at line $eof")

# Scale and write new file
open("geom.inp", "w") do f
    print(f, "*NODE, NSET=Nall\n")  # write the first line

    for i in 4:elem_def
        line = strip(lines[i])
        occursin("**", line) && break
        line_split = split(replace(line, "," => " "))
        data = to_float(line_split[2:end])
        scaled_coords = data #[scale, scale, zscale] .* data
        scaled_line = line_split[1] * join([", " * string(round(d,digits=6)) for d in scaled_coords])
        print(f, scaled_line * "\n")
    end

    print(f, "*ELEMENT, type=C3D4, ELSET=FLAP" * "\n")  # write element definition line

    for i in elem_def+1:eof
        if !occursin("*ELEMENT, type=C3D4,", lines[i])
            print(f, lines[i] * "\n")
        end
    end
end

println("File written to geom.inp")
