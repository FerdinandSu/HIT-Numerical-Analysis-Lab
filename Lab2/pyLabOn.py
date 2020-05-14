import pypandoc

class ReportContent(object):
    def __init__(self):
        super().__init__()







class ReportPlainContent(ReportContent):
    def __init__(self, content):
        super().__init__()
        self.content = content

    def __str__(self):
        return self.content + '\n'


class ReportTable(ReportContent):
    def __init__(self, headers, data, align=0):
        """
        :param headers: headers of the Table
        :param data: data of the table
        :param align align setting of the table:Left=-1,Right=1,Central(default)=0;
        :type align: list or int
        """
        super().__init__()
        self.Align = align if isinstance(align, list) else [align] * len(headers)
        self.Data = data
        self.Headers = headers
        self.__TableAlign = [' :----: |', ' ----:', ':-----|']

    def __str__(self):
        r = '|'
        width = len(self.Headers)
        for header in self.Headers:
            r += header + "|"

        r += '\n|'
        for align in self.Align:
            r += self.__TableAlign[align]
        r += '\n'
        for line in self.Data:
            r += '|'
            for cell in line:
                r += str(cell) + '|'
            r += '\n'

        return r


class ReportParagraph(ReportContent):

    def __init__(self, level, title):
        """
        Initialize a report paragraph with Title-level and Title
        :param level: Title-level
        :param title: Title
        :type level: int
        """
        super().__init__()
        self.SubContents = []

        self.Level = level
        self.Title = title

    def __str__(self):
        r = '#' * self.Level + ' ' + self.Title + '\n\n'
        for subContent in self.SubContents:
            r += str(subContent) + '\n'

        return r

    def add_sub_paragraph(self, title):
        sub_paragraph = ReportParagraph(self.Level + 1, title)
        self.SubContents.append(sub_paragraph)
        return sub_paragraph

    def add_table(self,headers, data, align=0):
        table=ReportTable(headers,data,align)
        self.SubContents.append(table)
        return table

    def add_plain_content(self,content):
        plain_content=ReportPlainContent(content)
        self.SubContents.append(plain_content)
        return plain_content

class Report(ReportParagraph):
    def __init__(self,title,path=''):
        super().__init__(1,title)
        self.MarkdownFilePath=path+'/'+title+'.md'
        self.PDFPath = path + '/' + title + '.pdf'
        self.MarkdownFile=open(self.MarkdownFilePath,"w+")

    def save_to_file(self):
        file = open(self.MarkdownFilePath, "w+")
        file.write(str(self))
        file.close()

    def compile(self,engine='xelatex',font="SimSun"):
        pypandoc.convert_file(self.MarkdownFilePath, 'pdf', outputfile=self.PDFPath,
                              extra_args=['--pdf-engine='+engine, '-V', 'mainfont='+font])