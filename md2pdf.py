""" Script to convert markdown to pdf """

from markdown_pdf import MarkdownPdf, Section

# Create a MarkdownPdf object
pdf = MarkdownPdf()

# Load Markdown content from a file or string
with open("input.md", "r") as f:
  markdown_content = f.read()

# Add a section to the PDF
pdf.add_section(Section(markdown_content))

# Save the PDF
pdf.save("output.pdf")
