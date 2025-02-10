import csv
import math
import time
import pandas as pd

def parse_tsv(file_path):
    return pd.read_csv(file_path, sep='\t', skiprows=1)

def calculate_distance(ra1, dec1, ra2, dec2):
    ra1, dec1, ra2, dec2 = map(math.radians, [ra1, dec1, ra2, dec2])
    delta_ra = ra1 - ra2
    delta_dec = dec1 - dec2
    a = math.sin(delta_dec / 2) ** 2 + math.cos(dec1) * math.cos(dec2) * math.sin(delta_ra / 2) ** 2
    return 2 * math.asin(math.sqrt(a))

def adjust_ra(ra, star_ra):
    new_ra = star_ra - ra
    if new_ra >= 360: new_ra -= 360
    elif new_ra < 0: new_ra += 360
    return new_ra

def adjust_dec(dec, star_dec):
    new_dec = star_dec - dec
    if new_dec > 90: new_dec -= 180
    elif new_dec < -90: new_dec += 180
    return new_dec

def filter_stars(stars, ra, dec, fov_h, fov_v):
    filtered_stars = []
    for _, star in stars.iterrows():
        star_ra = star["ra_ep2000"]
        star_dec = star["dec_ep2000"]
        range_ra = abs(adjust_ra(ra, star_ra))
        range_dec = abs(adjust_dec(dec, star_dec))
        if range_ra <= fov_h / 2 and range_dec <= fov_v / 2:
            filtered_stars.append(star)
    return filtered_stars

def sort_stars_by_brightness(stars):
    n = len(stars)
    for i in range(n):
        for j in range(0, n - i - 1):
            if stars[j]["brightness"] > stars[j + 1]["brightness"]:
                # Swap if the current star is dimmer than the next star
                stars[j], stars[j + 1] = stars[j + 1], stars[j]
    return stars

def top_stars(stars, n):
    return stars[:n]

def write_to_csv(stars, ra, dec):
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    output_file = f"{timestamp}.csv"
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["ID", "RA", "DEC", "Brightness", "AngularDifference"])
        for star in stars:
            new_ra = adjust_ra(ra, star["ra_ep2000"])
            new_dec = adjust_dec(dec, star["dec_ep2000"])
            writer.writerow([star["source_id"], new_ra, new_dec, star["brightness"], star["distance_rad"]])

def main():
    # Input parameters
    tsv_file = "cleaned_stars.tsv"
    ra = float(input("Enter RA: "))
    dec = float(input("Enter DEC: "))
    fov_h = float(input("Enter FOV Horizontal (fov_h): "))
    fov_v = float(input("Enter FOV Vertical (fov_v): "))
    n = int(input("Enter the number of stars (N): "))

    # Parse the TSV file
    stars = parse_tsv(tsv_file)

    # phot_g_mean_mag appears to show the magnitude/apparent brightness of the star.
    stars.rename(columns={"phot_g_mean_mag": "brightness"}, inplace=True)  # Replace with the actual brightness column name from the TSV file
    stars["distance_rad"] = [calculate_distance(ra, dec, star_ra, star_dec) for star_ra, star_dec in zip(stars["ra_ep2000"], stars["dec_ep2000"])]

    # Filter stars within the field of view
    filtered_stars = filter_stars(stars, ra, dec, fov_h, fov_v)

    # Sort stars by distance
    sorted_stars = sort_stars_by_brightness(filtered_stars)

    # Select the top N the brightest stars
    top_n_stars = sorted_stars[:n]

    # Write the result to a CSV file
    write_to_csv(top_n_stars, ra, dec)

if __name__ == "__main__":
    main()
