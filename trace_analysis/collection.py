import os
import sys
import numpy as np
import re
import contextlib
import joblib.parallel
import multiprocessing
from tqdm import tqdm
from collections import UserList
# from collections.abc import MutableSequence

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


# Following code snippet from https://stackoverflow.com/questions/24983493/tracking-progress-of-joblib-parallel-execution/58936697#58936697
@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar given as argument"""
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()

def fun2(fun, obj, *args, **kwargs):
    fun_result = fun(obj, *args, **kwargs)
    # if hasattr(obj, __getstate__):
    #
    #     return (obj.__getstate__, fun_result)
    # if hasattr(__dict__)
    # else:
    return (obj, fun_result)

class Collection(UserList):
    def __init__(self, data=None, use_parallel_processing=True, number_of_cores=None, parallel_processing_kwargs=dict(verbose=0)):
        # self.data = data

        if number_of_cores is None:
            number_of_cores = multiprocessing.cpu_count()

        if data is None:
            super(Collection, self).__setattr__('data', [])
        else:
            super(Collection, self).__setattr__('data', data)
        # super(Collection, self).__init__(data)

        # super(Collection, self).__setattr__('parallel', 4)
        super(Collection, self).__setattr__('use_parallel_processing', use_parallel_processing)
        super(Collection, self).__setattr__('number_of_cores', number_of_cores)
        super(Collection, self).__setattr__('parallel_processing_kwargs', parallel_processing_kwargs)

        # data_types = {type(datum) for datum in self.data}
        # data_types.discard(type(None))
        # if len(data_types) > 1:
        #     raise NotImplementedError('At the moment only a single data type within a Collection is allowed')
        # super(Collection, self).__setattr__('data_type', list(data_types)[0])

    @property
    def dict_without_data(self):
        d = self.__dict__.copy()
        d.pop('data')
        # d.pop('data_type')
        return d

    @property
    def parallel(self):
        self.use_parallel_processing = True
        return self

    @property
    def serial(self):
        self.use_parallel_processing = False
        return self

    def map(self, fun):
        if not self.use_parallel_processing or len(self.data) == 1:
            print('Serial processing')

            def f(*args, **kwargs):
                with HiddenPrints():
                    # output = [getattr(datum, item)(*args, **kwargs) if datum is not None else None
                    #           for datum in tqdm(self.data, position=0, leave=True)]
                    output = [fun(datum, *args, **kwargs) if datum is not None else None
                              for datum in tqdm(self.data, position=0, leave=True)]
                for o in output:
                    if o is not None:
                        return Collection(output, **self.dict_without_data)

            return f
        else:
            print('Parallel processing')

            def f(*args, **kwargs):
                with HiddenPrints():
                    # , require='sharedmem')
                    # with tqdm_joblib(tqdm(self.data, position=0, leave=True)):
                    #     output = joblib.parallel.Parallel(self.number_of_cores, **self.parallel_processing_kwargs)\
                    #         (joblib.parallel.delayed(getattr(datum, item))(*args, **kwargs) for datum in self.data)

                    with tqdm_joblib(tqdm(self.data, position=0, leave=True)):
                        results = joblib.parallel.Parallel(self.number_of_cores, **self.parallel_processing_kwargs) \
                            (joblib.parallel.delayed(fun2)(fun, datum, *args, **kwargs) for datum in self.data)
                    #
                    #
                    # # with Pool() as pool:
                    # #     output = pool.starmap(fun2, [(fun, datum, args, kwargs) for datum in self.data])
                    # #
                    #
                    # output = Parallel(self.number_of_cores, verbose=10)\
                    #     (delayed(fun)(datum, *args, **kwargs) for datum in self.data)
                # output = Parallel(self.number_of_cores)(
                #     delayed(getattr(File, item))(datum, *args, **kwargs) if datum is not None else None for datum in
                #     tqdm(self.data))
                output = []
                for datum, (obj, out) in zip(self.data, results):
                    if hasattr(obj, '__dict__'):
                        obj.__dict__ = {key: value for key, value in obj.__dict__.items() if
                                        not hasattr(value, '_do_not_update')}
                        if hasattr(datum, '__setstate__'):
                            datum.__setstate__(obj.__getstate__())
                        elif hasattr(datum, '__dict__'):
                            datum.__dict__.update(obj.__dict__)
                    output.append(out)

                for o in output:
                    if o is not None:
                        return Collection(output, **self.dict_without_data)

            return f

    def __getattr__(self, item):
        if 'data' not in self.__dict__.keys():
            raise AttributeError
        # if item == 'data':
        #     super(Collection, self).__getattribute__(item)
        # if inspect.ismethod(getattr(self.files[0], item))

        # first_not_none = None
        for datum in self.data:
            if datum is not None:
                first_not_none = datum
                break

        # if first_not_none is not None:
        #     raise ValueError('Empty collection')

        if callable(getattr(first_not_none, item)):
            fun = getattr(type(self.data[0]), item)
            return self.map(fun)
        else:
            return Collection([getattr(datum, item) if datum is not None else None for datum in self.data], **self.dict_without_data)

    def __setattr__(self, key, value):
        # if key == 'data':
        #     super(Collection, self).__setattr__(key, value)
        if key in self.__dict__.keys():
            super(Collection, self).__setattr__(key, value)
        # elif hasattr(self.data[0], key):
        else:
            if len(self.data) == 0:
                raise IndexError('Collection is empty')
            for datum in self.data:
                if datum is not None:
                    setattr(datum, key, value)

    def set_values(self, key, values):
        if len(self.data) == 0:
            raise IndexError('Collection is empty')
        for datum, value in zip(self.data, values):
            setattr(datum, key, value)

    def __getitem__(self, item):
        # type(self) is used instead of Collection to enable proper returning in child classes
        if isinstance(item, int):
            return self.data[item]
        elif isinstance(item, np.ndarray) or isinstance(item, list) or isinstance(item, Collection):
            return type(self)(np.array(self.data)[item].tolist(), **self.dict_without_data)
            # if isinstance(data, list) and len(data) == 1:
            #     return data[0]
            # else:
            #     return Collection(data)
        else: # For slicing
            return type(self)(self.data[item], **self.dict_without_data)


    def __delitem__(self, key):
        self.data.pop(key)

    def __len__(self):
        return len(self.data)

    def __setitem__(self, key, value):
        self.data[key] = value

    def insert(self, index, object):
        self.data.insert(index, object)

    @property
    def str(self):
        return Collection([str(datum) for datum in self.data], **self.dict_without_data)

    def regex(self, pattern):
        p = re.compile(pattern.replace('\\', r'\\'))
        return [True if p.search(string) else False for string in self.data]

    def not_none(self):
        return Collection([datum for datum in self.data if datum is not None], **self.dict_without_data)
