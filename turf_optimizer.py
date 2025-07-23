import math
import numpy as np
import json
from shapely.geometry import Polygon, box
from shapely.affinity import rotate, translate
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
from io import BytesIO
import tkinter as tk
from tkinter import filedialog
import os
import sys
import traceback

print("==== SCRIPT STARTED ====")

def optimize_roll_rotation(coords, roll_width):
    """
    Brute-force search over angles 0…359° (no offset).
    Return a dictionary for the layout with minimum waste area.
    """
    original_poly = Polygon(coords)
    centroid = original_poly.centroid
    best_layout = None
    min_waste = float("inf")
    for angle in range(360):
        rotated_poly = rotate(original_poly, angle, origin=centroid)
        # Use zero offset coverage for simplicity in this optimization
        placed, used, strips = _coverage_with_offset(rotated_poly, roll_width, offset=0.0)
        if not strips:
            continue
        waste = placed - used
        if waste < min_waste:
            min_waste = waste
            best_layout = {
                "rotated_polygon": rotated_poly,
                "used_area": used,
                "placed_area": placed,
                "waste": waste,
                "roll_polygons": strips,
                "angle": angle,
                "original_area": rotated_poly.area,
                "num_rolls": len(strips)
            }
    return best_layout

def pareto_frontier(coords, roll_width):
    """
    Compute Pareto-optimal solutions (angle with minimal waste for each roll count).
    Returns a dict {roll_count: {...angle, waste, ...}, ...} for non-dominated entries.
    """
    poly = Polygon(coords)
    centroid = poly.centroid
    candidates = []
    for ang in range(360):
        rot = rotate(poly, ang, origin=centroid)
        placed, used, strips = _coverage_with_offset(rot, roll_width, offset=0.0)
        if not strips:
            continue
        candidates.append({
            "angle": ang,
            "used_area": used,
            "placed_area": placed,
            "waste": placed - used,
            "num_rolls": len(strips)
        })
    # Best (lowest waste) for each roll count
    best_per_roll = {}
    for res in candidates:
        r = res["num_rolls"]
        if r not in best_per_roll or res["waste"] < best_per_roll[r]["waste"]:
            best_per_roll[r] = res
    # Filter out dominated solutions across roll counts
    pareto_solutions = {}
    solutions = list(best_per_roll.values())
    for sol in solutions:
        dominated = False
        for other in solutions:
            if (other["num_rolls"] <= sol["num_rolls"] and other["waste"] <= sol["waste"] 
                    and (other["num_rolls"] < sol["num_rolls"] or other["waste"] < sol["waste"])):
                dominated = True
                break
        if not dominated:
            pareto_solutions[sol["num_rolls"]] = sol
    return pareto_solutions

def generate_turf_layout(data: dict, angle_step: int = 1, offset_steps: int = 10) -> bytes:
    """
    Build a PDF (bytes) with turf roll layout solutions. 
    If multiple Pareto-optimal solutions exist, a summary page with a Pareto frontier chart 
    and one page per solution is included. If only one solution is optimal, only a single layout page is produced.
    
    Parameters
    ----------
    data : dict 
        Expected keys:
          "coords": [(x0,y0), (x1,y1), ...]  (polygon vertices in metres) 
          OR "vertices": [{"lat": ..., "lng": ...}, ...] (polygon vertices in lat/lon, WGS84)
          "roll_width": <roll width in metres> OR "grass_width_m": <roll width in metres>
          Optional cost parameters:
          "grass_type": <name of grass/turf type>
          "grass_price_isk_per_m2": <cost of grass per square meter>
          "extras_detail": [ { "name": ..., "unit": ..., "price": ... }, ... ] (optional extra items)
          "sand_rate_kg_per_m2": <kilograms of sand needed per square meter>
          "sand_total_kg": <total kg of sand>
          "sand_total_cost_isk": <total cost of sand>
    Returns
    -------
    bytes
        PDF file content.
    """
    # Debug print statements
    print("Input data:", json.dumps(data, indent=2))
    print("Extras detail:", json.dumps(data.get("extras_detail", []), indent=2))
    
    # 1. Ensure coordinates are in metres (project lat/lon to local plane if needed)
    if "coords" not in data:
        if "vertices" not in data:
            raise ValueError("Need either 'coords' (metres) or 'vertices' (lat/lon).")
        verts = data["vertices"]
        if len(verts) < 3:
            raise ValueError("Polygon must have at least 3 vertices.")
        # Compute approximate local projection using centroid lat/lon
        lat0 = sum(v["lat"] for v in verts) / len(verts)
        lon0 = sum(v["lng"] for v in verts) / len(verts)
        deg2rad = math.pi / 180.0
        cos_lat0 = math.cos(lat0 * deg2rad)
        R = 6_371_000.0  # Earth radius in meters
        data["coords"] = [
            (
                (v["lng"] - lon0) * deg2rad * cos_lat0 * R,  # x in m
                (v["lat"] - lat0) * deg2rad * R             # y in m
            )
            for v in verts
        ]
    coords = data["coords"]
    if len(coords) < 3:
        raise ValueError("coords must contain at least 3 vertices")
    
    # Normalize orientation - align longest side with x-axis
    # Find the longest side of the polygon
    longest_length = 0
    longest_idx = 0
    for i in range(len(coords)):
        j = (i + 1) % len(coords)  # Next vertex (wrapping around to first)
        dx = coords[j][0] - coords[i][0]
        dy = coords[j][1] - coords[i][1]
        length = math.sqrt(dx*dx + dy*dy)
        
        if length > longest_length:
            longest_length = length
            longest_idx = i
    
    # Calculate angle of longest side with x-axis
    i = longest_idx
    j = (longest_idx + 1) % len(coords)
    dx = coords[j][0] - coords[i][0]
    dy = coords[j][1] - coords[i][1]
    angle_rad = math.atan2(dy, dx)
    angle_deg = math.degrees(angle_rad)
    
    # Create a polygon and rotate it to align longest side with x-axis
    original_polygon = Polygon(coords)
    normalized_polygon = rotate(original_polygon, -angle_deg, origin='centroid')
    
    # Get normalized coordinates and create polygon
    normalized_coords = list(normalized_polygon.exterior.coords)[:-1]  # Remove last point which duplicates first
    
    print("Original longest side angle:", angle_deg)
    print("After normalization, longest side should be aligned with x-axis")
    
    polygon = normalized_polygon
    if not polygon.is_valid:
        polygon = polygon.buffer(0)  # fix self-intersections if any
    if polygon.area <= 0:
        raise ValueError("Invalid polygon geometry (area must be positive)")

    shape_area = polygon.area
    # Use provided roll width (in metres). Accept legacy key "grass_width_m" as alias.
    roll_width = float(data["roll_width"] if "roll_width" in data else data.get("grass_width_m", 0))
    data["roll_width"] = roll_width  # ensure it's recorded in data

    # Extract cost data with appropriate key mappings
    grass_type = data.get("grass_type")
    grass_cost_per_m2 = data.get("grass_price_isk_per_m2")
    
    # Extract extras from the extras_detail list
    extras_detail = data.get("extras_detail", [])
    
    # Debug print for extras_detail
    print("Received extras_detail:", extras_detail)
    
    # Find glue info in extras
    glue_item = next((item for item in extras_detail if item.get("name") == "Lim"), None)
    glue_cost_per_m = glue_item.get("price") if glue_item else None
    
    # Handle sand info
    sand_kg_per_m2 = data.get("sand_rate_kg_per_m2")
    sand_total_kg = data.get("sand_total_kg")
    sand_cost_total = data.get("sand_total_cost_isk")
    sand_cost_per_kg = None
    if sand_total_kg and sand_cost_total:
        sand_cost_per_kg = sand_cost_total / sand_total_kg if sand_total_kg > 0 else 0

    # Debug print for sand and glue parameters
    print("Sand parameters:", {
        "sand_kg_per_m2": sand_kg_per_m2,
        "sand_total_kg": sand_total_kg,
        "sand_cost_total": sand_cost_total,
        "sand_cost_per_kg": sand_cost_per_kg
    })
    print("Glue parameters:", {
        "glue_cost_per_m": glue_cost_per_m
    })

    # 2. Search for layout solutions by varying angle and offset
    best_per_roll = {}
    for ang in range(0, 180, angle_step):
        rot_poly = rotate(polygon, ang, origin="centroid")
        minx, miny, maxx, maxy = rot_poly.bounds
        for k in range(offset_steps):
            off = (k / offset_steps) * roll_width
            placed, used, strips = _coverage_with_offset(rot_poly, roll_width, off)
            if not strips:
                continue
            r = len(strips)
            waste = placed - used
            # Record the best (lowest waste) layout for this roll count
            if r not in best_per_roll or waste < best_per_roll[r]["waste"]:
                best_per_roll[r] = {
                    "angle": ang,
                    "offset": off,
                    "num_rolls": r,
                    "waste": waste,
                    "used": used,
                    "placed": placed,
                    "strips": strips    # store strip polygons for later glue calculation
                }

    if not best_per_roll:
        raise RuntimeError("No valid turf layout found for the given polygon and roll width.")

    # 3. Filter out dominated solutions (truly Pareto-optimal entries)
    layouts = sorted(best_per_roll.values(), key=lambda x: x["num_rolls"])
    pareto_layouts = []
    for layout in layouts:
        rolls_i = layout["num_rolls"]
        waste_i = layout["waste"]
        dominated = False
        for other in layouts:
            if (other["num_rolls"] <= rolls_i and other["waste"] <= waste_i 
                    and (other["num_rolls"] < rolls_i or other["waste"] < waste_i)):
                dominated = True
                break
        if not dominated:
            pareto_layouts.append(layout)
    # Sort Pareto-optimal layouts by roll count for consistent ordering
    pareto_layouts.sort(key=lambda x: x["num_rolls"])

    # 4. Determine the single best layout (minimum waste) for summary/labeling
    best_layout = min(pareto_layouts, key=lambda x: x["waste"])
    # Prepare a record for the best layout with all needed fields
    best_record = {
        "angle": best_layout["angle"],
        "offset": best_layout.get("offset", 0.0),
        "original_area": shape_area,
        "used_area": best_layout["used"],
        "placed_area": best_layout["placed"],
        "waste": best_layout["waste"],
        "num_rolls": best_layout["num_rolls"]
    }
    # Calculate glue length for the best layout
    glue_len_best = calculate_glue_length(best_layout["strips"])

    # 5. Generate PDF output
    pdf_buffer = BytesIO()
    with PdfPages(pdf_buffer) as pdf:
        if len(pareto_layouts) > 1:
            # Multiple solutions: include summary page with Pareto table & chart
            pareto_rows = []
            for sol in pareto_layouts:
                angle = sol["angle"]
                waste = sol["waste"]
                offset = sol["offset"]
                # Recalculate glue for each Pareto-optimal solution
                glue_len = calculate_glue_length(sol["strips"])
                rolls = sol["num_rolls"]
                pareto_rows.append((rolls, angle, waste, glue_len))
            # Build summary page
            _page_summary(pdf, best_record, glue_len_best, pareto_rows)
            # One page per Pareto-optimal layout
            for i, (rolls, angle, waste, glue_len) in enumerate(pareto_rows):
                # Get the offset from the original layout object
                offset = pareto_layouts[i]["offset"]
                _page_layout(pdf, polygon, roll_width, rolls, angle, offset, waste, glue_len, 
                             best_record, grass_type, grass_cost_per_m2, glue_cost_per_m, 
                             sand_kg_per_m2, sand_cost_per_kg, extras_detail, 
                             sand_total_kg, sand_total_cost_isk=sand_cost_total)
        else:
            # Only one optimal solution: just output the single layout page
            sol = pareto_layouts[0]
            # Compute glue length for this solution
            glue_len = calculate_glue_length(sol["strips"])
            _page_layout(pdf, polygon, roll_width, sol["num_rolls"], sol["angle"], sol["offset"], sol["waste"], glue_len, 
                         best_record, grass_type, grass_cost_per_m2, glue_cost_per_m, 
                         sand_kg_per_m2, sand_cost_per_kg, extras_detail, 
                         sand_total_kg, sand_total_cost_isk=sand_cost_total)
    pdf_bytes = pdf_buffer.getvalue()
    pdf_buffer.close()
    return pdf_bytes

def calculate_glue_length(roll_polygons) -> float:
    """
    Calculate the total length of shared edges between adjacent strips (i.e., glue seams in metres).
    """
    glue = 0.0
    for i, p1 in enumerate(roll_polygons):
        for j, p2 in enumerate(roll_polygons):
            if i >= j:
                continue
            inter = p1.intersection(p2)
            if inter.is_empty:
                continue
            if inter.geom_type == "LineString":
                glue += inter.length
            elif inter.geom_type == "MultiLineString":
                glue += sum(line.length for line in inter.geoms)
    return glue

def _coverage_with_offset(rot_poly: Polygon, roll_w: float, offset: float):
    """
    Compute turf coverage for a rotated polygon at a given vertical offset.
    Returns (placed_area, used_area, strips) where:
      - placed_area: total area of turf strips placed (including waste overlap),
      - used_area: area of turf that covers the polygon (shape area covered),
      - strips: list of strip polygons covering the shape.
    """
    minx, miny, maxx, maxy = rot_poly.bounds
    # Generate horizontal strip positions from bottom (miny) upwards, with a starting offset
    y_positions = np.arange(miny - offset, maxy, roll_w)
    placed = 0.0
    used = 0.0
    strips = []
    for y in y_positions:
        strip_rect = box(minx - 10, y, maxx + 10, y + roll_w)  # slightly wider than shape to fully cover
        inside = strip_rect.intersection(rot_poly)
        if inside.is_empty:
            continue
        # Area of shape covered by this strip:
        used += inside.area
        
        # Check if the intersection is a single polygon or multiple disconnected polygons
        if inside.geom_type == 'Polygon':
            # Process single polygon as before
            sx_min, _, sx_max, _ = inside.bounds
            placed_strip = box(sx_min, y, sx_max, y + roll_w)
            placed += placed_strip.area
            strips.append(placed_strip)
        else:
            # Handle MultiPolygon or GeometryCollection by processing each component separately
            geometries = []
            if inside.geom_type == 'MultiPolygon':
                geometries = list(inside.geoms)
            elif inside.geom_type == 'GeometryCollection':
                geometries = [geom for geom in inside.geoms if geom.geom_type == 'Polygon']
            
            for geom in geometries:
                # Create a strip for each separate polygon
                sx_min, _, sx_max, _ = geom.bounds
                placed_strip = box(sx_min, y, sx_max, y + roll_w)
                placed += placed_strip.area
                strips.append(placed_strip)
                
    return placed, used, strips

def _page_summary(pdf, best: dict, glue_len_opt: float, pareto_rows: list):
    """
    Creates the first PDF page: a summary with optimal layout info, Pareto table, and Pareto frontier chart.
    """
    fig = plt.figure(figsize=(8.27, 11.69))  # A4 size
    fig.patch.set_facecolor("white")
    # Remove summary text block
    
    # Pareto table (top-right)
    tbl_ax = fig.add_axes([0.08, 0.80, 0.84, 0.17])
    tbl_ax.axis("off")
    col_labels = ["Rolls", "Angle (°)", "Waste (m²)", "Glue (m)"]
    body = [[r, f"{a:.0f}", f"{w:.2f}", f"{g:.2f}"] for (r, a, w, g) in pareto_rows]
    table = tbl_ax.table(cellText=body, colLabels=col_labels, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.3)
    # Pareto chart (bottom)
    chart_ax = fig.add_axes([0.08, 0.10, 0.84, 0.55])
    chart_ax.grid(True, alpha=0.3)
    # Plot waste vs rolls
    chart_ax.plot([r for r, *_ in pareto_rows], 
                  [w for *_, w, __ in pareto_rows], marker="o")
    chart_ax.set_xlabel("Number of rolls")
    chart_ax.set_ylabel("Waste area (m²)")
    chart_ax.set_title("Pareto Frontier")
    # Save summary page
    pdf.savefig(fig)
    plt.close(fig)

def _page_layout(pdf,
                 polygon: Polygon,
                 roll_width: float,
                 rolls: int,
                 angle: float,
                 offset: float,
                 waste: float,
                 glue_len: float,
                 best_record: dict,
                 grass_type: str,
                 grass_cost_per_m2,
                 glue_cost_per_m,
                 sand_kg_per_m2,
                 sand_cost_per_kg,
                 extras_detail: list,
                 sand_total_kg=None,
                 sand_total_cost_isk=None):
    """
    Create a layout page for a given Pareto-optimal solution.
    Includes the layout diagram (left) and a detailed info block and cost table (right).
    """
    # For consistent naming, use sand_cost_total consistently inside this function
    sand_cost_total = sand_total_cost_isk
    
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.patch.set_facecolor("white")
    # Main layout plot on left half
    ax = fig.add_axes([0.08, 0.25, 0.60, 0.60])  # [left, bottom, width, height]
    
    # Rotate shape to the specified angle and translate to (0,0)
    rotated_poly = rotate(polygon, angle, origin="centroid")
    minx, miny, maxx, maxy = rotated_poly.bounds
    rotated_poly_orig = rotated_poly  # Store original rotated polygon before translation
    rotated_poly = translate(rotated_poly, xoff=-minx, yoff=-miny)
    width = maxx - minx
    height = maxy - miny
    
    # Draw accurate turf roll strips from coverage calculation using provided offset
    _, _, strip_polygons = _coverage_with_offset(rotated_poly_orig, roll_width, offset)
    
    for strip in strip_polygons:
        strip_trans = translate(strip, xoff=-minx, yoff=-miny)
        xs, ys = strip_trans.exterior.xy
        ax.fill(xs, ys, facecolor='lightgray', edgecolor='green', linewidth=1.0, alpha=0.6)
    
    # Draw shape area on top of strips (blue fill with green border)
    x_coords, y_coords = rotated_poly.exterior.xy
    ax.fill(x_coords, y_coords, color="blue", alpha=0.4, edgecolor="green")
    
    # Annotate side lengths (in meters) on the shape
    verts = list(rotated_poly.exterior.coords)
    for i in range(len(verts) - 1):
        x1, y1 = verts[i]
        x2, y2 = verts[i + 1]
        midx, midy = (x1 + x2) / 2, (y1 + y2) / 2
        length = math.hypot(x2 - x1, y2 - y1)
        ax.text(midx, midy, f"{length:.2f} m",
                fontsize=9, color="black", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))
                
    # Label the area in the center of the shape
    ax.text(rotated_poly.representative_point().x, rotated_poly.representative_point().y,
            f"{rotated_poly.area:.2f} m²", fontsize=12, color="white",
            ha="center", va="center",
            bbox=dict(boxstyle="round,pad=0.3", fc="black", ec="black"))
            
    # Format plot
    ax.set_facecolor("white")
    ax.set_aspect("equal", "box")
    ax.set_xlim(0, width + roll_width * 0.1)
    ax.set_ylim(0, height + roll_width * 0.1)
    ax.axis('off')  # Remove axis and ticks

    # Header text at top-center
    header_line1 = f"Rolls used: {rolls} (each {roll_width:.2f} m wide)"
    if grass_type:
        header_line1 += f" – {grass_type}"
    fig.text(0.5, 0.96, header_line1, ha="center", va="top", fontsize=12)
    fig.text(0.5, 0.92, f"Optimal angle: {angle:.0f}°", ha="center", va="top", fontsize=12)
    # Right-side info block (Area, Waste, Total)
    total_area = polygon.area + waste  # area + waste
    info_lines = [
        f"Area: {polygon.area:.2f} m²",
        f"Waste: {waste:.2f} m²",
        f"Total: {total_area:.2f} m²"
    ]
    fig.text(0.72, 0.85, "\n".join(info_lines), ha="left", va="top", fontsize=11)

    # Cost analysis title and table at bottom
    has_cost_data = (grass_cost_per_m2 or glue_cost_per_m or (sand_kg_per_m2 and sand_cost_per_kg) or extras_detail)
    print("\n---- COST TABLE GENERATION ----")
    print("has_cost_data:", has_cost_data)
    print("grass_cost_per_m2:", grass_cost_per_m2)
    print("glue_cost_per_m:", glue_cost_per_m)
    print("sand_kg_per_m2:", sand_kg_per_m2)
    print("sand_cost_per_kg:", sand_cost_per_kg)
    print("extras_detail:", json.dumps(extras_detail, indent=2))
    
    if has_cost_data:
        fig.text(0.5, 0.18, "Cost Analysis", ha="center", va="bottom", fontsize=12, fontweight="bold")
        cost_ax = fig.add_axes([0.08, 0.05, 0.84, 0.12])
        cost_ax.axis("off")
        # Build cost table rows
        rows = []
        total_cost_val = 0.0
        
        # Grass/turf cost
        if grass_cost_per_m2:
            grass_area = polygon.area + waste  # total turf area used (m²)
            grass_cost = grass_area * grass_cost_per_m2
            total_cost_val += grass_cost
            rows.append([
                grass_type or "Grass",
                f"{grass_area:.2f} m²",
                f"{grass_cost_per_m2:,.0f} kr",
                f"{grass_cost:,.0f} kr"
            ])
            
        # Find sand price in extras_detail if not provided directly
        if sand_kg_per_m2 and not sand_cost_per_kg:
            for extra in extras_detail:
                if extra.get("name") == "Sandur":
                    sand_cost_per_kg = extra.get("price")
                    print("Found sand price in extras:", sand_cost_per_kg)
                    break
                    
        # Process extras from extras_detail list
        processed_items = set()
        
        print("\n-- Processing Sand Items --")
        # First handle sand specifically
        for extra in extras_detail:
            name = extra.get("name", "")
            print(f"Checking if '{name}' is Sandur")
            if name == "Sandur":
                unit = extra.get("unit", "")
                # Use the values directly from the input data if available
                print("Found Sandur item:", extra)
                print("sand_total_kg:", sand_total_kg)
                print("sand_cost_total:", sand_cost_total)
                
                if sand_total_kg and sand_cost_total:
                    rows.append([
                        name,
                        f"{sand_total_kg:,.2f} {unit}",
                        f"{sand_cost_per_kg:,.0f} kr",
                        f"{sand_cost_total:,.0f} kr"
                    ])
                    total_cost_val += sand_cost_total
                    print("Added Sandur row with total values")
                elif sand_kg_per_m2 and sand_cost_per_kg:
                    # Calculate if specific values weren't provided
                    sand_amount = polygon.area * sand_kg_per_m2
                    sand_cost = sand_amount * sand_cost_per_kg
                    total_cost_val += sand_cost
                    rows.append([
                        name,
                        f"{sand_amount:,.2f} {unit}",
                        f"{sand_cost_per_kg:,.0f} kr",
                        f"{sand_cost:,.0f} kr"
                    ])
                    print("Added Sandur row with calculated values")
                else:
                    print("Not enough data to add Sandur row")
                processed_items.add("Sandur")
                break
        else:
            print("No Sandur item found in extras_detail")
        
        print("\n-- Processing Other Items --")
        # Then handle all other items
        for extra in extras_detail:
            name = extra.get("name", "")
            print(f"Processing extra: {name}")
            if name in processed_items or name == grass_type or name == "Nafn" or name == "Name":
                print(f"Skipping {name} (already processed or special item)")
                continue
                
            unit = extra.get("unit", "")
            price = extra.get("price")
            print(f"{name}: unit={unit}, price={price}")
            
            # Handle glue/tape items
            if name in ["Lim", "Limbordi"] and price is not None:
                # Use glue length for any glue-related items
                cost_per_m = price
                if name == "Lim" and glue_cost_per_m is not None:
                    cost_per_m = glue_cost_per_m
                    
                total_cost = glue_len * cost_per_m
                total_cost_val += total_cost
                rows.append([
                    name,
                    f"{glue_len:.2f} {unit}",
                    f"{cost_per_m:,.0f} kr",
                    f"{total_cost:,.0f} kr"
                ])
                print(f"Added {name} row with cost {total_cost:.0f} kr")
                processed_items.add(name)
            elif price is not None:  # Other extras with price info
                # Skip items with zero price and no other information
                if price == 0 and not (unit or "amount" in extra):
                    print(f"Skipping {name} (zero price)")
                    continue
                    
                amount = extra.get("amount", 1)  # Default to 1 if no amount specified
                total_cost = amount * price
                total_cost_val += total_cost
                rows.append([
                    name,
                    f"{amount} {unit}",
                    f"{price:,.0f} kr",
                    f"{total_cost:,.0f} kr"
                ])
                print(f"Added generic item {name} with cost {total_cost:.0f} kr")
                processed_items.add(name)
            
        print("\n-- Final Rows --")
        print(rows)
        
        # Final total row
        if total_cost_val and rows:
            rows.append(["Total Cost", "", "", f"{total_cost_val:,.0f} kr"])
        
        # Create table
        col_labels = ["Description", "Amount", "Unit Price", "Total Price"]
        table = cost_ax.table(cellText=rows, colLabels=col_labels, loc="center", cellLoc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.2)
        # Set table grid lines and bold styling for total row
        for (row, col), cell in table.get_celld().items():
            cell.set_edgecolor("black")
            cell.set_linewidth(0.4)
            text = cell.get_text().get_text()
            # Bold the "Total Cost" row (identify by label)
            if text.strip().lower() == "total cost" or text.strip().lower() == f"{total_cost_val:,.0f} kr".lower():
                cell.get_text().set_fontweight("bold")
    # Save layout page
    pdf.savefig(fig)
    plt.close(fig)

def main():
    """
    Load and process all JSON files in the examples folder, generating a PDF for each.
    """
    import glob
    from datetime import datetime
    try:
        print("\n\n===== MAIN FUNCTION STARTED =====")
        print("Current working directory:", os.getcwd())

        # 1) Decide which files to process
        if len(sys.argv) > 1:  # python turf_optimizer.py path/to/file.json
            # Robust path handling for CLI input
            cli_path = sys.argv[1]
            event_path = os.path.normpath(cli_path)
            tried_paths = [event_path]
            if not os.path.exists(event_path):
                # Try with .json extension
                if not event_path.lower().endswith('.json'):
                    event_path_json = event_path + '.json'
                    tried_paths.append(event_path_json)
                    if os.path.exists(event_path_json):
                        event_path = event_path_json
                # Try with .json.json extension
                if not os.path.exists(event_path):
                    if not event_path.lower().endswith('.json.json'):
                        event_path_jsonjson = event_path + '.json'
                        tried_paths.append(event_path_jsonjson)
                        if os.path.exists(event_path_jsonjson):
                            event_path = event_path_jsonjson
            if not os.path.exists(event_path):
                raise FileNotFoundError(f"Input file not found. Tried: {tried_paths}")
            input_files = [event_path]
            print(f"Using input file from CLI: {input_files[0]}")
        else:
            # Find all .json or .json.json files in examples/
            input_files = sorted(glob.glob(os.path.join("examples", "*.json")))
            input_files += sorted(glob.glob(os.path.join("examples", "*.json.json")))
            if not input_files:
                raise FileNotFoundError(
                    "No input files found in examples/. Run: python turf_optimizer.py examples/event_sample_rectangle.json"
                )
            print(f"Found {len(input_files)} input files in examples/: {input_files}")

        for event_path in input_files:
            print(f"\n--- Processing {event_path} ---")
            with open(event_path, "r") as f:
                raw = json.load(f)

            if isinstance(raw, dict) and "body" in raw and isinstance(raw["body"], str):
                data = json.loads(raw["body"])
            else:
                data = raw

            print("Generating turf layout...")
            pdf_bytes = generate_turf_layout(data, angle_step=1, offset_steps=10)

            # Output name: use input file name and timestamp
            base = os.path.splitext(os.path.basename(event_path))[0]
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_name = f"turf_layout_{base}_{timestamp}.pdf"

            print(f"Saving PDF to {output_name}...")
            with open(output_name, 'wb') as f:
                f.write(pdf_bytes)

            print(f"Success! Generated PDF saved to {output_name}")

        print("===== MAIN FUNCTION COMPLETED =====\n")

    except Exception as e:
        print("\n===== ERROR IN MAIN FUNCTION =====")
        print(f"Error type: {type(e).__name__}")
        print(f"Error message: {str(e)}")
        print("Traceback:")
        traceback.print_exc()
        print("===== END OF ERROR TRACEBACK =====\n")

    finally:
        print("Press Enter to exit...")
        input()

if __name__ == "__main__":
    print("Calling main() function...")
    main()
    print("Main function returned, script ending.")

print("==== SCRIPT ENDED ====")
