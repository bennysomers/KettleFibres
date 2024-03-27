column_boundaries = [3, 8, 15, 20, 28, 40, 51, 62,
                     73, 83, 94, 107, 118, 128, 138, 147, 158, 166, ]

# Function to parse the text and convert it to CSV format


def parse_text_to_csv(text):
    # Split the text into lines
    lines = text.strip().split('\n')

    # Define the headers for the CSV
    headers = [
        "Element", "A", "Z", "Iso", "Spin", "Maxw.(n-g)", "error", "(n-g)", "error", "(n-el)", "error", "(n-a)", "error", "(n-p)", "error", "(n-f)", "error", "Iso"
    ]

    # Initialize an empty list to store CSV rows
    csv_rows = [','.join(headers)]

    for i in range(10, len(lines)):
        last_index = 0
        line = [None] * len(column_boundaries)
        for j in range(0, len(column_boundaries)):
            line[j] = lines[i][last_index:column_boundaries[j]]
            last_index = column_boundaries[j]
            line[j] = line[j].strip()
        # Append the columns to the CSV row
        csv_row = ','.join(line)
        csv_rows.append(csv_row)

    # Join all CSV rows with newline characters
    csv_text = '\n'.join(csv_rows)

    return csv_text

# Function to read input from a file


def read_input_file(file_path):
    with open(file_path, 'r') as file:
        text = file.read()
    return text

# Main function


def main():
    # Path to the input file
    # Change this to the actual path of your input file
    input_file_path = 'capture_parameters/thermal_resonance_ranges.txt'

    # Read input from the file
    input_text = read_input_file(input_file_path)

    # Parse the input text and convert it to CSV format
    csv_text = parse_text_to_csv(input_text)

    # Print the CSV text
    print(csv_text)


if __name__ == "__main__":
    main()
