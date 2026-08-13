from PySide2.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QPlainTextEdit, QLabel
)
from PySide2.QtCore import Qt
import sys
import io


from papylio.gui.common_layouts import HelpDialog


class ScriptBox(QWidget):

    def __init__(self, parent=None):
        super().__init__(parent)

        self.parent = parent
        self._file = None

        # ------------------------------------------------------------
        # Script editor
        # ------------------------------------------------------------

        self.script_edit = QPlainTextEdit()
        self.script_edit.setPlaceholderText(
            "Paste Python/Papylio code here..."
        )

        # ------------------------------------------------------------
        # Output window
        # ------------------------------------------------------------

        self.output = QPlainTextEdit()
        self.output.setReadOnly(True)

        # ------------------------------------------------------------
        # Buttons
        # ------------------------------------------------------------

        self.run_button = QPushButton("Run")
        self.clear_button = QPushButton("Clear")
        self.help_button = QPushButton("Help")

        self.run_button.clicked.connect(self.run_script)
        self.clear_button.clicked.connect(self.clear_output)
        self.help_button.clicked.connect(self.show_script_help)

        button_layout = QHBoxLayout()
        button_layout.addWidget(self.run_button)
        button_layout.addWidget(self.clear_button)
        button_layout.addStretch()
        button_layout.addWidget(self.help_button)

        # ------------------------------------------------------------
        # Main layout
        # ------------------------------------------------------------

        layout = QVBoxLayout(self)

        layout.addWidget(QLabel("Python script:"))
        layout.addWidget(self.script_edit)

        layout.addLayout(button_layout)

        layout.addWidget(QLabel("Output:"))
        layout.addWidget(self.output)

        self.setDisabled(True)


    # ================================================================
    # Current file
    # ================================================================

    @property
    def file(self):
        return self._file

    @file.setter
    def file(self, file):
        self._file = file

        if file is None:
            self.setDisabled(True)
        else:
            self.setDisabled(False)


    # ================================================================
    # Run script
    # ================================================================

    def run_script(self):

        code = self.script_edit.toPlainText()

        if not code.strip():
            return

        self.output.clear()
        self.output.appendPlainText("Running script...\n")

        # Objects available to the script
        namespace = {
            "self": self.parent,
            "file": self.file,
            "experiment": self.parent.experiment,
        }

        # Capture print() output
        stdout = io.StringIO()
        stderr = io.StringIO()

        old_stdout = sys.stdout
        old_stderr = sys.stderr

        sys.stdout = stdout
        sys.stderr = stderr

        try:
            exec(code, namespace)

        except Exception as e:
            # Put the error into our captured output
            print(f"{type(e).__name__}: {e}")

        finally:
            # VERY IMPORTANT: restore normal Python output
            sys.stdout = old_stdout
            sys.stderr = old_stderr

        # Show captured output in the GUI
        output_text = stdout.getvalue()
        error_text = stderr.getvalue()

        if output_text:
            self.output.appendPlainText(output_text)

        if error_text:
            self.output.appendPlainText(error_text)

        self.output.appendPlainText("\nScript finished.")

    # ================================================================
    # Clear output
    # ================================================================

    def clear_output(self):
        self.output.clear()


    # ================================================================
    # Help
    # ================================================================

    def show_script_help(self):

        help_text = """
        <html>
        <body style="font-family: sans-serif; font-size: 10pt;">

        <h2>Custom Scripts</h2>

        <p>
        Paste and execute Python/Papylio scripts.
        </p>

        <p>
        The following objects are available:
        </p>

        <ul>
          <li><b>file</b> - currently selected file</li>
          <li><b>experiment</b> - current experiment</li>
          <li><b>self</b> - the MainWindow</li>
        </ul>

        <h3>Example</h3>

        <pre>
        #show current mapping
        mapping_file = experiment.files[0]
        figure, axis = mapping_file.show_image()
        mapping_file.mapping.show(axis=axis, show_source=True)
        figure.show()
        </pre>
        
        <h3>Note<\h3>
        
        <p>
         Custom code may overwrite earlier GUI pipeline results.
         It may lead to mismatches and or errors in the stored data.
         The user is assumed to have a proper insight in Papylio script.
        </p>
        
        

        </body>
        </html>
        """

        self.help_dialog = HelpDialog(self, help_text)
        self.help_dialog.show()