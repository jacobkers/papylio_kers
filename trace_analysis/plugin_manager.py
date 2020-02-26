import importlib
import pkgutil

class PluginManager:
    def __init__(self):
        self.plugins_module = importlib.import_module('trace_analysis.plugins')
        self.plugin_names = [pluginname for _, pluginname, ispkg in pkgutil.walk_packages(self.plugins_module.__path__) if ispkg]
        self.plugins = [importlib.import_module('trace_analysis.plugins.'+plugin_name) for plugin_name in self.plugin_names]

    def get_class_plugins(self, class_name):
        return tuple([getattr(plugin, class_name) for plugin in self.plugins if hasattr(plugin, class_name)])



