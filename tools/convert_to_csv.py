import csv

# Open the text file and read its contents.
with open('Offsets_0.txt', 'r') as infile:
    content = infile.read()

# Split the content into individual number strings (assuming whitespace separation).
numbers = content.split()

# Open a new CSV file for writing.
with open('Offsets_0.csv', 'w', newline='') as outfile:
    writer = csv.writer(outfile)
    
    # Write the header row.
    writer.writerow(['Index', 'Value'])
    
    # Write each number with its corresponding index.
    for index, number in enumerate(numbers):
        writer.writerow([index, number])

print("CSV file has been created successfully as 'Offsets_0.csv'.")
