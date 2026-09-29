# Shared look + helpers for every FigureN/compose.jl (include at the top). Provides:
#   LEGEND                       legend entries (labels describe the biology, not PhysiCell type names)
#   zoom(xlim, ylim)             Panel transform: crop a snapshot to a µm window, with a scale bar
#   prefix_ids(tag)              Panel transform for pdftocairo schematics (avoids id clashes)
#   stack_rows(rows)             stack rows of different heights into one SVG (Montage grids are uniform)
#   snapshot_minutes(file), snapshot_cells(file)   read a snapshot's time (min) and cell count
using Montage

# Legend entries are written out (not read from PhysiCell's legend.svg) so labels describe the
# biology: the proliferating neoplastic cells use PhysiCell's "CAF" cell type.
LEGEND = [
    ("Epithelial",                 "grey"),
    ("Neoplastic (proliferating)", "yellow"),
    ("Basement membrane",          "red"),
]

# PhysiCell draws snapshots at 1 px = 1 µm: svg x = x - X_MIN, svg y = SVG_H - (y - Y_MIN)
# (Duct_Dev_Figs domain is ±400 µm; SVG_H includes PhysiCell's 56 px caption band)
const X_MIN, Y_MIN, SVG_H = -400, -400, 856

"""
    zoom(xlim, ylim; bar = 50)

Panel transform that crops a PhysiCell snapshot to the window `xlim` × `ylim` (µm), with a
thin frame and a `bar` µm scale bar. The crop drops PhysiCell's caption, so put the time in the title.
"""
function zoom(xlim, ylim; bar = 50)
    w, h = xlim[2] - xlim[1], ylim[2] - ylim[1]
    sx, sy = xlim[1] - X_MIN, SVG_H - (ylim[2] - Y_MIN)        # window's top-left corner in svg px
    return function (svg)
        body = replace(svg, r"^.*?<svg\b.*?>"s => "", r"</svg>\s*$" => "")   # strip root tag
        return """
        <svg width="$w" height="$h" xmlns="http://www.w3.org/2000/svg">
        <svg width="$w" height="$h" overflow="hidden"><g transform="translate($(-sx),$(-sy))">$body</g></svg>
        <rect x="$(w - bar - 8)" y="$(h - 10)" width="$bar" height="3" fill="black"/>
        <text x="$(w - bar / 2 - 8)" y="$(h - 14)" text-anchor="middle" font-family="Arial" font-size="9">$bar µm</text>
        <rect x="0" y="0" width="$w" height="$h" fill="none" stroke="black" stroke-width="1"/>
        </svg>"""
    end
end

# Panel transform for schematics converted with pdftocairo: every such file reuses the same ids
# (glyph-0-1, clip-1, …), so prefix them with `tag` before several share one figure. Also turns
# `xlink:href` into plain `href`, since nesting drops the root tag that declared the xlink namespace.
prefix_ids(tag) = svg -> replace(svg, "id=\"" => "id=\"$tag-", "xlink:href=\"#" => "href=\"#$tag-", "url(#" => "url(#$tag-")

"""
    stack_rows(rows; width = 1260, pad = 12)

Arrange SVGs (file paths or SVG text) in `rows`: items in a row share one height and together span
`width`; rows stack top to bottom at their natural heights (unlike Montage's uniform grid cells).
Returns the SVG text.
"""
function stack_rows(rows; width = 1260, pad = 12)
    io = IOBuffer()
    y = pad
    for row in rows
        texts = [Montage._svgSource(item) for item in row]
        dims = Montage._svgDimensions.(texts)
        h = (width - (length(row) + 1) * pad) / sum(w / hh for (w, hh) in dims)   # row height that fills the width
        x = pad
        for (text, (w, hh)) in zip(texts, dims)
            println(io, Montage._nestedSVG(text, x, y, w * h / hh, h))
            x += w * h / hh + pad
        end
        y += h + pad
    end
    return """
    <svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$y" viewBox="0 0 $width $y">
    <rect x="0" y="0" width="$width" height="$y" fill="white"/>
    $(String(take!(io)))</svg>"""
end

# Simulated time of a snapshot (minutes), read from PhysiCell's "Current time" caption
function snapshot_minutes(file)
    m = match(r"Current time: (\d+) days, (\d+) hours, and ([\d.]+) minutes", read(file, String))
    return 1440 * parse(Int, m[1]) + 60 * parse(Int, m[2]) + parse(Float64, m[3])
end

# Cell count of a snapshot, read from PhysiCell's "N agents" caption
snapshot_cells(file) = parse(Int, match(r"(\d+) agents", read(file, String))[1])
