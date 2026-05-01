"""
    gml_parser.jl — Parse GML files from the Internet Topology Zoo dataset.
    
    Extracts nodes (with id, label, latitude, longitude, internal flag)
    and edges (source, target) from .gml formatted files.
"""

struct GMLNode
    id::Int
    label::String
    latitude::Float64
    longitude::Float64
    internal::Bool
end

struct GMLEdge
    source::Int
    target::Int
end

struct GMLGraph
    name::String
    nodes::Vector{GMLNode}
    edges::Vector{GMLEdge}
end

function parse_gml(filepath::String)
    lines = readlines(filepath)

    nodes = GMLNode[]
    edges = GMLEdge[]
    graph_name = basename(filepath)
    graph_name = replace(graph_name, ".gml" => "")

    i = 1
    while i <= length(lines)
        line = strip(lines[i])

        if startswith(line, "Network ")
            m = match(r"Network\s+\"(.+)\"", line)
            if m !== nothing
                graph_name = m.captures[1]
            end
        end

        if line == "node ["
            node_id = -1
            label = ""
            lat = NaN
            lon = NaN
            internal = true
            i += 1
            depth = 1
            while i <= length(lines) && depth > 0
                nline = strip(lines[i])
                if endswith(nline, "[")
                    depth += 1
                elseif nline == "]"
                    depth -= 1
                    if depth == 0
                        break
                    end
                end
                if depth == 1
                    if startswith(nline, "id ")
                        node_id = parse(Int, split(nline)[2])
                    elseif startswith(nline, "label ")
                        m = match(r"label\s+\"(.*)\"", nline)
                        if m !== nothing
                            label = m.captures[1]
                        else
                            label = split(nline, " ", limit=2)[2]
                        end
                    elseif startswith(nline, "Latitude ")
                        lat = parse(Float64, split(nline)[2])
                    elseif startswith(nline, "Longitude ")
                        lon = parse(Float64, split(nline)[2])
                    elseif startswith(nline, "Internal ")
                        internal = parse(Int, split(nline)[2]) == 1
                    end
                end
                i += 1
            end
            if node_id >= 0
                push!(nodes, GMLNode(node_id, label, lat, lon, internal))
            end
        end

        if line == "edge ["
            src = -1
            tgt = -1
            i += 1
            depth = 1
            while i <= length(lines) && depth > 0
                eline = strip(lines[i])
                if endswith(eline, "[")
                    depth += 1
                elseif eline == "]"
                    depth -= 1
                    if depth == 0
                        break
                    end
                end
                if depth == 1
                    if startswith(eline, "source ")
                        src = parse(Int, split(eline)[2])
                    elseif startswith(eline, "target ")
                        tgt = parse(Int, split(eline)[2])
                    end
                end
                i += 1
            end
            if src >= 0 && tgt >= 0
                push!(edges, GMLEdge(src, tgt))
            end
        end

        i += 1
    end

    return GMLGraph(graph_name, nodes, edges)
end

function list_gml_files(data_dir::String)
    files = filter(f -> endswith(lowercase(f), ".gml"), readdir(data_dir))
    return sort([joinpath(data_dir, f) for f in files])
end
