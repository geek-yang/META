#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the meta module.
"""
import pytest

from meta import meta


def test_without_test_object():
    assert False


class TestMeta(object):
    @pytest.fixture
    def return_a_test_object(self):
        pass

    def test_meta(self, meta):
        assert False

    def test_with_error(self, meta):
        with pytest.raises(ValueError):
            pass
