# HTML template functions for GenerateReport.jl
function report_css()
    return """
:root{--bg:#f4f7f6;--card:#fff;--txt:#2c3e50;--accent:#3498db;--axis:#34495e;--border:#eaedf1;--side:280px}
html{scroll-behavior:smooth}
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,sans-serif;background:var(--bg);color:var(--txt);margin:0;display:flex;overflow-x:hidden}
.sidebar{width:var(--side);height:100vh;position:fixed;top:0;left:0;background:#fff;border-right:1px solid #e0e6ed;padding:30px 20px;box-sizing:border-box;overflow-y:auto;box-shadow:2px 0 10px rgba(0,0,0,.03);z-index:100;transition:transform .3s}
.sidebar h3{margin-top:0;color:#2980b9;border-bottom:2px solid #ecf0f1;padding-bottom:12px;font-size:16px;text-transform:uppercase;letter-spacing:1px;padding-right:20px}
.sidebar ul{list-style:none;padding:0;margin:0}.sidebar li{margin-bottom:4px}
.sidebar a{text-decoration:none;color:#576574;font-weight:600;font-size:13px;display:block;padding:8px 12px;border-radius:6px;transition:all .2s}
.sidebar a:hover{background:#f0f4f8;color:var(--accent);padding-left:16px}
.sidebar .section-header{font-size:11px;text-transform:uppercase;color:#95a5a6;letter-spacing:1.5px;margin:16px 0 6px 4px;font-weight:700}
.close-btn{position:absolute;top:25px;right:15px;background:none;border:none;font-size:24px;color:#95a5a6;cursor:pointer}.close-btn:hover{color:#e74c3c}
.open-btn{position:fixed;top:20px;left:20px;z-index:90;background:var(--card);border:1px solid var(--border);padding:10px 15px;border-radius:8px;cursor:pointer;box-shadow:0 2px 8px rgba(0,0,0,.05);font-weight:600;color:var(--txt);display:none;align-items:center;gap:8px}.open-btn:hover{background:#f8f9fa}
body.sidebar-closed .sidebar{transform:translateX(-100%)}body.sidebar-closed .open-btn{display:flex}body.sidebar-closed .main{margin-left:0;width:100%}
.main{margin-left:var(--side);padding:40px 60px;width:calc(100% - var(--side));box-sizing:border-box;transition:margin-left .3s,width .3s}
h1{text-align:center;color:var(--txt);margin-bottom:10px;font-weight:800;letter-spacing:-.5px}
.ts{text-align:center;color:#7f8c8d;font-size:14px;margin-bottom:30px;font-style:italic}
.bb{background:#e8f4f8;border-left:5px solid var(--accent);padding:20px;border-radius:8px;margin-bottom:40px;box-shadow:0 2px 4px rgba(0,0,0,.03)}
.bb h4{margin:0 0 15px;color:#2980b9;font-size:15px;text-transform:uppercase;letter-spacing:1px}
.bb ul{list-style:none;padding:0;margin:0;display:flex;flex-wrap:wrap;gap:15px}
.bb li{font-size:13px;background:#fff;padding:8px 15px;border-radius:5px;border:1px solid #dcebf0}
.filter-info{background:#fff3cd;border-left:5px solid #ffc107;padding:15px;border-radius:8px;margin-bottom:40px}
.filter-info h4{margin:0 0 8px;color:#856404;font-size:14px}
.filter-info p{margin:4px 0;font-size:13px;color:#856404}
h2{border-bottom:3px solid var(--accent);padding-bottom:10px;margin-top:60px;color:var(--txt);scroll-margin-top:40px;text-transform:capitalize}
.exp-section{background:linear-gradient(135deg,#667eea 0%,#764ba2 100%);color:#fff;padding:15px 25px;border-radius:10px;margin-top:50px;margin-bottom:10px;font-size:20px;font-weight:700;letter-spacing:.5px;scroll-margin-top:30px}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(240px,1fr));gap:20px;margin-bottom:40px}
.card{background:var(--card);border:1px solid var(--border);border-radius:12px;padding:15px;text-align:center;box-shadow:0 4px 6px rgba(0,0,0,.04);transition:transform .2s,box-shadow .2s}
.card:hover{transform:translateY(-3px);box-shadow:0 8px 15px rgba(0,0,0,.08)}
.card img{max-width:100%;height:auto;display:block;margin:0 auto;border-radius:6px}
.card .label{margin-top:12px;font-weight:600;font-size:14px;color:var(--txt);background:#f8f9fa;padding:6px;border-radius:6px;border:1px solid #e2e8f0}
.card .cc{font-size:11px;color:#7f8c8d;margin-top:4px}
.card .metrics{font-size:11px;color:#576574;margin-top:6px;background:#f0f4f8;padding:6px 8px;border-radius:4px;text-align:left;line-height:1.6}
.card .metrics span{display:inline-block;margin-right:8px}
.tc{background:var(--card);border-radius:12px;box-shadow:0 4px 15px rgba(0,0,0,.04);padding:20px;margin-bottom:50px;overflow-x:auto}
.tt{width:100%;border-collapse:separate;border-spacing:0}
.tt th,.tt td{border:1px solid var(--border);padding:12px;text-align:center;vertical-align:middle}
.corner{border:none!important;background:transparent!important}
.ax{background:var(--axis);color:#fff;font-weight:700;letter-spacing:1px;text-transform:uppercase;font-size:14px}
.xax{padding:12px;border-top-left-radius:8px;border-top-right-radius:8px}
.yax{width:40px;border-top-left-radius:8px;border-bottom-left-radius:8px}
.yax span{writing-mode:vertical-rl;transform:rotate(180deg);display:inline-block;white-space:nowrap}
.ch,.rh{background:#f8f9fa;color:var(--txt);font-weight:700}.rh{width:80px}
.tt img{max-width:180px;height:auto;display:block;margin:0 auto;border-radius:4px;transition:transform .2s}
.tt img:hover{transform:scale(1.08);box-shadow:0 4px 12px rgba(0,0,0,.15);z-index:10;position:relative}
.tt .cell-metrics{font-size:10px;color:#576574;margin-top:4px;line-height:1.5}
.empty{color:#a0aec0;font-size:13px;font-style:italic;background:#fafbfc}
@page{size:landscape;margin:10mm}
@media print{.sidebar,.open-btn{display:none!important}.main{margin-left:0!important;padding:0!important;width:100%!important}.grid{grid-template-columns:repeat(5,1fr)!important;gap:10px!important}.card img{max-width:140px!important}.tt img{max-width:110px!important}h2{page-break-after:avoid}}
"""
end
function report_js()
    return """<script>function toggleSidebar(){document.body.classList.toggle("sidebar-closed")}</script>"""
end
"""Format BM metrics as an HTML snippet for display in cards/table cells."""
function format_metrics_html(metrics::Dict{String,Float64}; compact::Bool=false)
    parts = String[]
    haskey(metrics, "IC") && push!(parts, "IC=$(round(metrics["IC"]; digits=3))")
    haskey(metrics, "area_frac_change") && push!(parts, "ΔA=$(round(100*metrics["area_frac_change"]; digits=1))%")
    haskey(metrics, "max_indent_depth") && push!(parts, "Depth=$(round(metrics["max_indent_depth"]; digits=3))")
    haskey(metrics, "indent_extent") && push!(parts, "Ext=$(round(100*metrics["indent_extent"]; digits=1))%")
    haskey(metrics, "lobe_amp") && push!(parts, "Lobe=$(round(metrics["lobe_amp"]; digits=3))")
    haskey(metrics, "roughness") && push!(parts, "Rough=$(round(metrics["roughness"]; digits=3))")
    haskey(metrics, "max_displacement") && push!(parts, "MaxD=$(round(metrics["max_displacement"]; digits=1))μm")
    isempty(parts) && return ""
    if compact
        return """<div class="cell-metrics">$(join(parts, " | "))</div>"""
    else
        return """<div class="metrics">$(join(["<span>$p</span>" for p in parts]))</div>"""
    end
end
function build_card_html(svg_filename::String, label::String, cell_count::Int;
                         is_baseline::Bool=false, metrics::Dict{String,Float64}=Dict{String,Float64}())
    cls = is_baseline ? "card baseline-card" : "card"
    bl = is_baseline ? " (Baseline)" : ""
    metrics_html = format_metrics_html(metrics)
    return """<div class="$(cls)"><img src="$(svg_filename)" alt="$(label)"><div class="label">$(label)$(bl)</div><div class="cc">cells: $(cell_count)</div>$(metrics_html)</div>\n"""
end
function build_grid_html(cards::Vector{String})
    return "<div class=\"grid\">\n" * join(cards) * "</div>\n"
end
function build_2d_table_html(display_y::String, display_x::String,
                             all_y::Vector{String}, all_x::Vector{String},
                             matrix::Dict{Tuple{String,String},Tuple{String,Int,Dict{String,Float64}}})
    nc, nr = length(all_x), length(all_y)
    h = """<div class="tc"><table class="tt"><thead>\n"""
    h *= """<tr><th colspan="2" rowspan="2" class="corner"></th><th colspan="$(nc)" class="ax xax">$(display_x)</th></tr>\n<tr>\n"""
    for xv in all_x
        h *= """<th class="ch">$(xv)</th>\n"""
    end
    h *= """</tr></thead>\n<tbody>\n"""

    for (y_idx, yv) in enumerate(all_y)
        h *= "<tr>\n"
        if y_idx == 1
            h *= """<td class="yax" rowspan="$(nr)"><span>$(display_y)</span></td>\n"""
        end
        h *= """<th class="rh">$(yv)</th>\n"""
        for xv in all_x
            if haskey(matrix, (yv, xv))
                fname, cc, metrics = matrix[(yv, xv)]
                metrics_html = format_metrics_html(metrics; compact=true)
                h *= """<td><img src="$(fname)" alt="$(display_y)_$(yv)_$(display_x)_$(xv)"><div class="cc">cells: $(cc)</div>$(metrics_html)</td>\n"""
            else
                h *= """<td class="empty">N/A</td>\n"""
            end
        end
        h *= "</tr>\n"
    end
    h *= "</tbody></table></div>\n"
    return h
end