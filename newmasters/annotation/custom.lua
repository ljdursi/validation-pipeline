function mean(vals)
    local sum=0
    for i=1,#vals do
        sum = sum + vals[i]
    end
    return sum / #vals
end

function sum(vals)
    local s=0
    for i=1,#vals do
        s = s + vals[i]
    end
    return s
end

function div2(a, b)
    if(a == 0) then return "0.0" end
    return string.format("%.9f", (a + 0) / b)
end

function strandedratio(a, b)
    return string.format("%.9f", sum(a) / b)
end

function split(str, sep)
    local sep, fields = sep or ":", {}
    local pattern = string.format("([^%s]+)", sep)
    str:gsub(pattern, function(c) fields[#fields+1] = c end)
    return fields
end

function by_priority(priority_list, set)
    for i=1,#priority_list do
        if set[priority_list[i]] then
            return priority_list[i]
        end
    end
    return "."
end

function prioritize_gencode(vals)
    local set = {}
    local v = split(vals, ",") or {}

    for i=1,#v do set[v[i]] = true end

    local priority_list = {"start_codon", "stop_codon", "gene", "transcript", "UTR", "CDS", "exon"}
    return by_priority(priority_list, set)
end
