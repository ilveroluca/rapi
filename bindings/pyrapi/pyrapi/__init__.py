
###############################################################################
# Copyright (c) 2014-2016 Center for Advanced Studies,
#                         Research and Development in Sardinia (CRS4)
# 
# Licensed under the terms of the MIT License (see LICENSE file included with the
# project).
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
###############################################################################

import pyrapi.rapi

def load_aligner(name="rapi_bwa"):
    """
    This function is used to load a specific aligner plugin.
    """
    return pyrapi.rapi
