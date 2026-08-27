#========================================================

# Figures for the SMC vs Flip vs Merge-split ensemble comparison

# Images save to working directory

#========================================================

#Required packages

library(redist)
library(sf)
library(spdep)
library(igraph)

sf_use_s2(FALSE) # avoids S2 errors from st_centroid/st_make_valid on fl25

#=================
# Load state data
#=================

data(fl25, package = "redist")
shp = st_as_sf(fl25)
shp = st_make_valid(shp)

#==========================
# Build the adjacency graph
#==========================

nb = poly2nb(shp, queen = TRUE) # contiguity list
adj_mat = nb2mat(nb, style = "B", zero.policy = TRUE) # binary adjacency matrix
g = graph_from_adjacency_matrix(adj_mat, mode = "undirected") # igraph object from that matrix

# County centroids, used as the graph's node positions
coords = suppressWarnings(st_coordinates(st_centroid(shp)))

#==================================
# State to graph figure (2 panels)
#==================================

png("state_graph_figure.png", width = 2000, height = 1100, res = 150, bg = "transparent")
par(mfrow = c(1, 2), bg = NA, mar = c(0, 0, 2, 0)) # small top margin so title sits close

# Panel A: the state / precincts
plot(st_geometry(shp), col = "cornflowerblue", border = "gray30")
title(main = "(A) State", col.main = "black", line = 0.3)

# Panel B: the underlying adjacency graph, overlaid on faint boundaries
plot(st_geometry(shp), col = NA, border = "gray30")
title(main = "(B) Graph", col.main = "black", line = 0.3)
plot(g, layout = coords, add = TRUE, rescale = FALSE, # layout = coords places nodes at real centroids, not an auto layout
     vertex.size = 1, vertex.color = "black", vertex.label = NA,
     edge.color = "blue")

dev.off()

#=============
# Shared setup
#=============

V(g)$name = as.character(seq_len(vcount(g))) # tag every vertex with its original id, so it survives subgraphing below
edges_idx = as_edgelist(g, names = FALSE) # plain matrix of vertex-id pairs, one row per edge

tree_edges = function(tg) matrix(as.integer(as_edgelist(tg, names = TRUE)), ncol = 2) # pull an igraph tree back out as an original-id edge matrix

spanning_tree = function(vids) {
  sg = induced_subgraph(g, vids) # the subgraph covering just these precincts
  E(sg)$weight = runif(ecount(sg)) # random edge weights
  mst(sg) # ...so the minimum spanning tree is effectively a random spanning tree
}

# Panel width/height (px) matching shp's true aspect ratio, so the map
# fills the panel with no letterboxing above/below
panel_dim = function(base_w) {
  bb = st_bbox(shp)
  asp = as.numeric((bb["ymax"] - bb["ymin"]) / (bb["xmax"] - bb["xmin"]))
  c(w = base_w, h = round(base_w * asp))
}

cex_main = 1.15
title_pad = 34 # px reserved per panel row for the title text

#=============================
# Merge-Split Steps (6 panels)
#=============================

set.seed(2)
pal = c("cornflowerblue", "orchid", "darkorange", "seagreen3")

district = kmeans(coords, centers = 3)$cluster #cluster precincts into districts
pair = c(2, 3) # the two districts that will be merged and resplit
d_merged = district
d_merged[d_merged %in% pair] = 4 # relabel both merge-target districts as one group, for coloring
merged_vids = which(district %in% pair) # precinct indices belonging to either district in `pair`

within_edges = edges_idx[district[edges_idx[, 1]] == district[edges_idx[, 2]], , drop = FALSE] # edges with both ends in the same district
tree1 = tree_edges(spanning_tree(which(district == 1)))
tree2 = tree_edges(spanning_tree(which(district == 2)))
tree3 = tree_edges(spanning_tree(which(district == 3)))

tree_merged = spanning_tree(merged_vids) # one new random spanning tree over the merged district
cut_e = E(tree_merged)[1] # the tree edge chosen to cut, for illustration
cut = as.integer(ends(tree_merged, cut_e)) # that edge's two endpoint precinct ids
tree_cut = delete_edges(tree_merged, cut_e) # removing it splits the tree into two pieces
comp = components(tree_cut)$membership # which of the two pieces each vertex now falls in
ids = as.integer(V(tree_cut)$name) # original precinct ids, in the same order as "comp"

new_district = district
new_district[ids[comp == 1]] = pair[1] # first piece becomes one new district
new_district[ids[comp == 2]] = pair[2] # second piece becomes the other

draw = function(d, ttl, edges = NULL, cut_pts = NULL) {
  plot(st_geometry(shp), col = pal[d], border = "white")
  title(ttl, line = 0.2, cex.main = cex_main)
  if (!is.null(edges) && nrow(edges)) # draw the supplied edges as line segments
    segments(coords[edges[, 1], 1], coords[edges[, 1], 2],
             coords[edges[, 2], 1], coords[edges[, 2], 2], lwd = 1.2)
  if (!is.null(cut_pts)) # highlight the cut edge in red and bolded for clear visual
    segments(coords[cut_pts[1], 1], coords[cut_pts[1], 2],
             coords[cut_pts[2], 1], coords[cut_pts[2], 2], col = "red", lwd = 3)
  points(coords, pch = 19, cex = 0.5)
}

pd = panel_dim(750)
png("recom_steps.png", width = pd["w"] * 3, height = (pd["h"] + title_pad) * 2,
    res = 170, pointsize = 20, bg = "transparent")
par(mfrow = c(2, 3), bg = NA, mar = c(0, 0, 1.6, 0))

plot(st_geometry(shp), border = "gray50")
title("(A) Graph", line = 0.2, cex.main = cex_main)
segments(coords[edges_idx[, 1], 1], coords[edges_idx[, 1], 2],
         coords[edges_idx[, 2], 1], coords[edges_idx[, 2], 2], col = "blue")
points(coords, pch = 19, cex = 0.5)

draw(district, "(B) District Graphs", within_edges)
draw(district, "(C) District Trees", rbind(tree1, tree2, tree3))
draw(d_merged, "(D) Merge & Sample New Tree", tree_edges(tree_merged))
draw(d_merged, "(E) Find Edge to Cut", tree_edges(tree_merged), cut)
draw(new_district, "(F) Split into Two Trees", tree_edges(tree_cut))

dev.off()

#=====================================
# Flip Algorithm (Boundary Swap steps)
#=====================================

set.seed(1)
district2 = kmeans(coords, centers = 2)$cluster
district_colors = c("lightblue", "salmon")

from = edges_idx[, 1]
to = edges_idx[, 2]
same_district = district2[from] == district2[to] # TRUE where an edge's two ends are in the same district
bnd = unique(c(from[!same_district], to[!same_district])) # precincts touching a district boundary

draw_districts = function(d, ttl) {
  plot(st_geometry(shp), col = NA, border = "gray80")
  title(ttl, line = 0.2, cex.main = cex_main)
  keep = d[from] == d[to] # only draw edges that stay within one district
  segments(coords[from[keep], 1], coords[from[keep], 2],
           coords[to[keep], 1], coords[to[keep], 2], col = "gray30", lwd = 1.5)
  points(coords, pch = 21, cex = 2, bg = district_colors[d], col = "black")
}

bnd_edges = cbind(from, to)[!same_district, , drop = FALSE] # edge list restricted to cross-district ("boundary") edges
deg = degree(g) # how many neighbors each precinct has in the full graph
best_edge = which.min(deg[bnd_edges[, 1]]) # a low-degree boundary precinct
selected = bnd_edges[best_edge, 1] # the precinct proposed to flip
partner = bnd_edges[best_edge, 2] # its neighbor on the other side of the boundary
target_district = district2[partner] # the district ""selected" would join if the flip is accepted

pd = panel_dim(700)
png("flip_steps.png", width = pd["w"] * 2, height = (pd["h"] + title_pad) * 2,
    res = 150, bg = "transparent")
par(mfrow = c(2, 2), bg = NA, mar = c(0, 0, 1.6, 0))

draw_districts(district2, '"Turn on" edges')

draw_districts(district2, "Gather connected components on boundaries")
symbols(coords[bnd, 1], coords[bnd, 2], circles = rep(0.03, length(bnd)),
        add = TRUE, inches = FALSE, lty = 2, fg = "black")

draw_districts(district2, "Select components and propose swaps")
symbols(coords[selected, 1], coords[selected, 2], circles = 0.04,
        add = TRUE, inches = FALSE, lty = 2, fg = "black")
segments(coords[selected, 1], coords[selected, 2],
         coords[partner, 1], coords[partner, 2], col = "red", lwd = 2)

new_district2 = district2
new_district2[selected] = target_district # apply the flip: move "selected" into "target_district"
draw_districts(new_district2, "Accept or reject the proposal")

dev.off()

#=====================
# SMC steps (4 panel)
#=====================

set.seed(3)
pal3 = c("plum", "darkorange", "steelblue", "seagreen3")
remaining = seq_len(vcount(g)) # precincts not yet assigned to a finished district
district_smc = rep(0, vcount(g)) # 0 = still unassigned

draw_smc = function(ttl) {
  cols = c("gray80", pal3)[district_smc + 1] # +1 shifts 0 (unassigned) to the first color, gray
  plot(st_geometry(shp), col = cols, border = "white")
  title(ttl, line = 0.2, cex.main = cex_main)
}

pd = panel_dim(550)
png("smc_steps.png", width = pd["w"] * 4, height = pd["h"] + title_pad,
    res = 160, pointsize = 20, bg = "transparent")
par(mfrow = c(1, 4), bg = NA, mar = c(0, 0, 1.6, 0))

draw_smc("(A) Initial map")

for (k in 1:3) {
  t = spanning_tree(remaining) # random spanning tree over whatever precincts remain
  
  # Try every possible single-edge cut, keep the most balanced one
  # (avoids peeling off small number of precincts that's hard to see)
  cuts = lapply(seq_len(ecount(t)), function(e) delete_edges(t, E(t)[e])) # every tree with one edge removed
  small_sizes = sapply(cuts, function(t2) min(table(components(t2)$membership))) # size of the smaller piece, for each candidate cut
  t2 = cuts[[which.max(small_sizes)]] # the cut whose smaller piece is as large as possible
  
  comp = components(t2)$membership # which of the two pieces each vertex falls in
  ids2 = as.integer(V(t2)$name) # original precinct ids, in the same order as "comp"
  small = as.integer(names(which.min(table(comp)))) # the label of the smaller piece
  new_ids = ids2[comp == small] # precinct ids in that smaller piece (rounds new district)
  
  district_smc[new_ids] = k
  remaining = setdiff(remaining, new_ids) # remove the newly-assigned precincts from the pool
  if (k == 3) district_smc[remaining] = k + 1 # last piece becomes the final district
  
  draw_smc(paste0("(", LETTERS[k + 1], ") Iteration ", k))
}

dev.off()