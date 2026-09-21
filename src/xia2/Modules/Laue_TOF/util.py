from __future__ import annotations

import logging
import time

import xia2.Driver.timing

xia2_logger = logging.getLogger(__name__)


def report_timing(fn):
    def wrap_fn(*args, **kwargs):
        start_time = time.time()
        result = fn(*args, **kwargs)
        xia2_logger.debug("\nTiming report:")
        xia2_logger.debug("\n".join(xia2.Driver.timing.report()))
        duration = time.time() - start_time
        # write out the time taken in a human readable way
        xia2_logger.info(
            "Processing took %s", time.strftime("%Hh %Mm %Ss", time.gmtime(duration))
        )
        return result

    return wrap_fn
