""" Script to convert markdown to pdf """

from markdown_pdf import MarkdownPdf, Section

# Create a MarkdownPdf object
pdf = MarkdownPdf()

filePath = r'C:\Users\U8019277\OneDrive - UniSQ\BlastOne_project'+r'\\'
fileName = r'safe_operation'
# Load Markdown content from a file or string
with open(filePath+fileName+r'.md', "r") as f:
  markdown_content = f.read()

# Add a section to the PDF
pdf.add_section(Section(markdown_content))

# Save the PDF
pdf.save(filePath+fileName+r'.pdf')

print('+++ PDF conversion completed +++')
